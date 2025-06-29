import ete3_extensions, pandas as pd, numpy as np, click, gzip, json, collections
from sklearn.model_selection import ParameterGrid
from sklearn.ensemble import RandomForestClassifier
from multiprocess import Pool

@click.command()
@click.option('-n', '--nexus', help='nexus file from EPHI', required=True)
@click.option('-m', '--matrix', help='reconstructed ancestral states by ETOKI', required=True)
@click.option('-w', '--weight', help='weight file. default: None', show_default=True, default=None)
@click.option('-p', '--prefix', help='prefix for output model', required=True)
@click.option('-l', '--leaf_only', help='prefix for output model', show_default=True, default=False, is_flag=True)
def train(prefix, matrix, nexus, weight, leaf_only) :
    pool = Pool(10)
    weights = {}
    if weight :
        weight_list = pd.read_csv(weight, sep='\t', header=None)
        for s, w in weight_list.values :
            try :
                weights[s] = float(w)
            except :
                pass

    # read nexus - node traits
    tre = ete3_extensions.read_nexus(nexus)[0]
    nodes = {}
    for node in tre.get_descendants() :
        if node.name :
            trait = node.annotations['state']
            if trait == 'ND' or (leaf_only and not node.is_leaf()) :
                continue
            if trait not in weights :
                weights[trait] = 1.
            weight = weights[trait] if node.is_leaf() else 1.
            nodes[node.name] = [trait, weight]
    state_code = {w:i for i, w in enumerate(weights)}
    nodes = {n:[state_code[info[0]], int(info[1]*100000000+0.5)] for n, info in nodes.items()}
    # read matrix - snps
    snp_matrix = []
    sites = []
    with gzip.open(matrix, 'rt') as fin :
        for line in fin :
            if line.startswith('##') :
                continue
            headers = line.strip().split('\t')
            header_ids = np.array([i for i, h in enumerate(headers) if i<2 or h in nodes])
            break
        for chunk in pd.read_csv(fin, header=None, sep='\t', usecols=header_ids, chunksize = 2048) :
            print(chunk.values[0,:2])
            idx = ~np.any(np.vectorize(len)(chunk.values[:, 2:]) > 1, 1)
            if np.sum(idx) :
                sites.append(chunk.values[idx, :2])
                snp_matrix.append(np.vectorize({'A':0, 'C':1, 'G':2, 'T':3, '-':4, '.':4}.get)(chunk.values[idx, 2:]).astype(np.uint8))
        sites = np.vstack(sites)
        snp_matrix = np.vstack(snp_matrix).T
        node_names = np.array(headers)[header_ids][2:]

    # seperate nodes into training, validation, and testing sets
    scores = []
    models = collections.defaultdict(list)
    for test_ite in np.arange(5) :
        idx = np.random.permutation(np.arange(node_names.size))
        mix_idx, test_idx = idx[:int(node_names.size*0.8+0.5)], idx[int(node_names.size*0.8+0.5):]
        test_name, test_mat = node_names[test_idx], snp_matrix[test_idx]
        test_y = pd.DataFrame([nodes[n] for n in test_name]).values
        v_step = mix_idx.size/5.

        for cv_ite in np.arange(5) :
            s, e = int(cv_ite*v_step+0.5), int((cv_ite+1)*v_step+0.5)
            train_idx, validate_idx = np.concatenate([mix_idx[:s], mix_idx[e:]]), mix_idx[s:e]
            train_name, train_mat = node_names[train_idx], snp_matrix[train_idx]
            validate_name, validate_mat = node_names[validate_idx], snp_matrix[validate_idx]
            train_y = pd.DataFrame([nodes[n] for n in train_name]).values
            validate_y = pd.DataFrame([nodes[n] for n in validate_name]).values
            # for each node, get top 30 closest nodes and traits
            train_x = get_x(train_mat, train_mat, train_y, len(state_code), ignore_self=True, pool=pool)
            validate_x = get_x(validate_mat, train_mat, train_y, len(state_code), pool=pool)

            m = find_model(train_x, train_y, validate_x, validate_y, pool)
            for score, key, _ in m :
                models[key].append(score)
        # for key, scores in models.items() :
        #     models[key] = np.mean(scores)
        print()
        best_model = max(models.items(), key=lambda m: np.mean(m[1][-5:]))
        best_param = json.loads(best_model[0])
        tv_name, tv_mat = node_names[mix_idx], snp_matrix[mix_idx]
        tv_y = pd.DataFrame([nodes[n] for n in tv_name]).values
        tv_x = get_x(tv_mat, tv_mat, tv_y, len(state_code), ignore_self=True, pool=pool)
        test_x = get_x(test_mat, tv_mat, tv_y, len(state_code), pool=pool)
        
        rfc = RandomForestClassifier(**best_param)
        rfc.fit(tv_x, tv_y[:, 0], sample_weight=tv_y[:, 1])
        score = rfc.score(test_x, test_y[:, 0], sample_weight=test_y[:, 1])
        scores.append(score)
        print(f'TEST: {best_model} --- {score}')
    
    final_model = max(models.items(), key=lambda m: np.mean(m[1]))
    print(final_model)
    final_param = json.loads(final_model[0])
    y = pd.DataFrame([nodes[n] for n in node_names]).values
    x = get_x(snp_matrix, snp_matrix, y, len(state_code), ignore_self=True, pool=pool)
    rfc = RandomForestClassifier(**final_param)
    rfc.fit(x, y[:, 0], sample_weight=y[:, 1])
    final_score = rfc.score(x, y[:, 0], sample_weight=y[:, 1])
    data = dict(model=rfc, params=final_model, scores=scores, final_score=final_score)
    import pickle
    with open(f'{prefix}.model.pickle', 'wb') as fout :
        pickle.dump(data, fout)





def find_model(train_x, train_y, validate_x, validate_y, pool) :
    params = dict(n_estimators=[250, 200, 150, 100], \
                  max_depth=[25, 20, 15, 10], \
                  max_features=[0.5, 0.4, 0.3, 0.2, "sqrt"])
    models = []
    for param, rfc in pool.imap_unordered(ite_find_model, [[param, train_x, train_y] for param in ParameterGrid(params)]) :
        score = rfc.score(validate_x, validate_y[:, 0], sample_weight=validate_y[:, 1])
        print(f'    {param} - {score}')
        models.append([score, param, rfc])
    print("\n")
    return models


def ite_find_model(data) :
    param, train_x, train_y = data
    rfc = RandomForestClassifier(**param)
    rfc.fit(train_x, train_y[:, 0], sample_weight=train_y[:, 1])
    return json.dumps(param, sort_keys=True), rfc


def ite_get_x(data) :
    mat, train_mat, train_y, min_size, n_feat, ignore_self = data
    x = []
    pp = np.sum(train_mat <4, 1)
    for i,m in enumerate(mat) :
        # if i % 100 == 0 : print(i)
        presence = np.sum((m == train_mat) & (m <4), 1)/pp
        idx = np.argsort(-presence)[1:n_feat+1] if ignore_self else np.argsort(-presence)[:n_feat]
        weights = presence[idx]/np.max(presence[idx])
        xx = np.zeros([n_feat, min_size])
        xx[np.arange(n_feat), train_y[idx, 0]] = weights
        x.append(xx.ravel())
    return x


def get_x(mat, train_mat, train_y, min_size, n_feat=30, ignore_self=False, pool=None) :
    x = []
    for xx in pool.map(ite_get_x, [[mat[ids], train_mat, train_y, min_size, n_feat, ignore_self] for ids in np.array_split(np.arange(mat.shape[0]), 10)]) :
        x.extend(xx)
    return x


if __name__ == '__main__' :
    train()
