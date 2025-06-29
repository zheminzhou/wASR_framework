import ete3_extensions, pandas as pd, numpy as np, click, gzip, sys
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from multiprocess import Pool

@click.command()
@click.option('-r', '--reference', help='reference file [needed for prediction]', required=True)
@click.option('-n', '--nexus', help='nexus file from EPHI', required=True)
@click.option('-m', '--matrix', help='reconstructed ancestral states by ETOKI', required=True)
@click.option('-w', '--weight', help='weight file. default: None', show_default=True, default=None)
@click.option('-p', '--prefix', help='prefix for output model', required=True)
@click.option('-N', '--n_neighbor', help='number of closest neighbors to consider', default=1, type=int)
def train(prefix, reference, matrix, nexus, weight, n_neighbor) :
    pool = Pool(10)
    weights = {}
    if weight :
        weight_list = pd.read_csv(weight, sep='\t', header=None)
        for s, w in weight_list.values :
            try :
                weights[s] = min(float(w), 1.)
            except :
                pass

    # read nexus - node traits
    tre = ete3_extensions.read_nexus(nexus)[0]
    nodes = {}
    for node in tre.get_descendants() :
        if node.name :
            traits = dict(zip(node.annotations['traits.list'], node.annotations['traits.prop']))
            ptraits = dict(zip(node.up.annotations['traits.list'], node.up.annotations['traits.prop'])) \
                if node.up else {'ND':1.}

            trait = max(traits.items(), key=lambda t:t[1])[0]
            if trait not in weights :
                weights[trait] = 1.
            weight = weights[trait] if node.is_leaf() else 1.
            nodes[node.name] = [traits, ptraits, weight]
    states = {t for traits in nodes.values() for t in traits[0].keys()}
    state_code = {s:i for i, s in enumerate(sorted(states))}
    node2 = {}
    for n, info in nodes.items() :
        n_state = np.zeros(len(states))
        for s, p in info[0].items() :
            n_state[state_code[s]] = p

        p_state = np.zeros(len(states))
        for s, p in info[1].items() :
            p_state[state_code[s]] = p
        
        node2[n] = [n_state, p_state, int(info[2]*100000000+0.5)]
    nodes = node2

    # read matrix - snps
    snp_matrix = []
    sites = []
    with gzip.open(matrix, 'rt') as fin :
        for line in fin :
            if line.startswith('##') :
                continue
            headers = line.strip().split('\t')
            header_ids = np.array([i for i, h in enumerate(headers) if i < 2 or h in nodes])
            break
        for chunk in pd.read_csv(fin, header=None, sep='\t', usecols=header_ids, chunksize = 4096) :
            print(chunk.values[0,:2])
            # if chunk.values[0, 1] > 1000000 :
            #     break
            idx = ~np.any(np.vectorize(len)(chunk.values[:, 2:]) > 1, 1)
            if np.sum(idx) :
                sites.append(chunk.values[idx, :2])
                snp_matrix.append(np.vectorize({'A':0, 'C':1, 'G':2, 'T':3, '-':4, '.':4}.get)(chunk.values[idx, 2:]).astype(np.uint8))
        sites = np.vstack(sites)
        snp_matrix = np.vstack(snp_matrix).T
        node_names = np.array(headers)[header_ids][2:]

    # seperate nodes into training, validation, and testing sets
    scores = []
    
    rand_idx = np.random.permutation(np.arange(node_names.size))
    v_step = rand_idx.size/5.
    
    y_score_flat, y_n_flat, y_a_flat, weight_flat = [], [], [], []
    for cv_ite in np.arange(5) :
        s, e = int(cv_ite*v_step+0.5), int((cv_ite+1)*v_step+0.5)
        train_idx, validate_idx = np.concatenate([rand_idx[:s], rand_idx[e:]]), rand_idx[s:e]
        train_name, train_mat = node_names[train_idx], snp_matrix[train_idx]
        validate_name, validate_mat = node_names[validate_idx], snp_matrix[validate_idx]
        train_y = pd.DataFrame([nodes[n] for n in train_name]).values
        validate_y = pd.DataFrame([nodes[n] for n in validate_name]).values

        y_pred_proba = np.array(get_x(validate_mat, train_mat, train_y, n_feat=n_neighbor, pool=pool))

        pred = np.argmax(y_pred_proba, 1)
        y_node = np.argmax([y for y in validate_y[:, 0]], 1)
        y_anc = np.argmax([y for y in validate_y[:, 1]], 1)
        
        score1 = np.sum((pred == y_node) * validate_y[:, 2]) / np.sum(validate_y[:, 2])
        score2 = np.sum(((pred == y_node) | (pred == y_anc)) * validate_y[:, 2]) / np.sum(validate_y[:, 2])

        scores.append([score1, score2])
        print(f'TEST_{cv_ite}: --- {score1} --- {score2}')
        sys.stdout.flush()
    
        y_n_bin = np.zeros_like(y_pred_proba, dtype=int)
        y_n_bin[np.arange(y_node.size), y_node] = 1

        idx = y_pred_proba[np.arange(y_node.size), y_node] > y_pred_proba[np.arange(y_node.size), y_anc]
        y_anc[idx] = y_node[idx]
        y_a_bin = np.zeros_like(y_pred_proba, dtype=int)
        y_a_bin[np.arange(y_node.size), y_anc] = 1
    
        # Calculate micro-average ROC curve and ROC area
        y_score_flat.append(y_pred_proba.ravel())
        y_n_flat.append(y_n_bin.ravel())
        y_a_flat.append(y_a_bin.ravel())
        weight_flat.append(np.repeat(validate_y[:, 2], y_pred_proba.shape[1]))
    
    data = dict(reference=reference, matrix=matrix, nexus=nexus, weights=weights,
                state_code=state_code, scores=np.array(scores))
    print(data)
    
    import pickle
    with open(f'{prefix}.model.pickle', 'wb') as fout :
        pickle.dump(data, fout)


    y_score_flat = np.concatenate(y_score_flat)
    y_n_flat = np.concatenate(y_n_flat)
    y_a_flat = np.concatenate(y_a_flat)
    weight_flat = np.concatenate(weight_flat)
    
    plt.figure(figsize=(12, 8))

    fpr, tpr, _ = roc_curve(y_n_flat, y_score_flat, sample_weight=weight_flat)
    roc_auc = auc(fpr, tpr)
    
    plt.plot(fpr, tpr, lw=2, label=f'Micro-average 1 (AUC = {roc_auc:.2f})')

    fpr2, tpr2, _ = roc_curve(y_a_flat, y_score_flat, sample_weight=weight_flat)
    roc_auc2 = auc(fpr2, tpr2)
    
    plt.plot(fpr2, tpr2, lw=2, label=f'Micro-average 2 (AUC = {roc_auc2:.2f})')
    
    # Plot settings
    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for Phylogenetic Trait Prediction')
    plt.legend(loc="lower right")
    
    # Save the plot
    output_file = f"{prefix}.roc_curves.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"ROC curves saved to {output_file}")
    
    plot_file = f'{prefix}.roc_curves.list'
    with open(plot_file, 'wt') as fout :
        p = [-1, -1]
        for x, y in zip(fpr, tpr) :
            if x - p[0] > 1e-4 or y - p[1] > 1e-4 :
                p = [x, y]
                fout.write(f'Curr_geo\t{x}\t{y}\n')
        if min(p) < 1 :
            fout.write(f'Curr_geo\t1\t1\n')
        fout.write('\n')
        p = [-1, -1]
        for x, y in zip(fpr2, tpr2) :
            if x - p[0] > 1e-4 or y - p[1] > 1e-4 :
                p = [x, y]
                fout.write(f'Trans_geo\t{x}\t{y}\n')
        if min(p) < 1 :
            fout.write(f'Trans_geo\t1\t1\n')



def ite_get_x(data) :
    mat, train_mat, train_y, n_feat, ignore_self = data
    n_feat2 = n_feat + 40
    x = []
    pp = (train_mat < 4)
    for i, m in enumerate(mat) :
        presence = np.sum((m == train_mat) & (m < 4), 1)/np.sum(pp * (m < 4), 1)
        idx2 = np.argsort(-presence)[1:n_feat2+1] if ignore_self else np.argsort(-presence)[:n_feat2]
        w = np.array([min(ww, 100000000) for ww in train_y[idx2, 2]])*(presence[idx2] - min(presence[idx2])+1e-8)
        i2 = np.argsort(-w)[:n_feat]
        res = np.sum(np.array(train_y[idx2[i2], 0].tolist()).T*w[i2], 1)
        x.append(res/np.sum(res))
    return x

def get_x(mat, train_mat, train_y, n_feat, ignore_self=False, pool=None) :
    x = []
    for xx in pool.map(ite_get_x, [[mat[ids], train_mat, train_y, n_feat, ignore_self] for ids in np.array_split(np.arange(mat.shape[0]), 10)]) :
        x.extend(xx)
    return x


if __name__ == '__main__' :
    train()
