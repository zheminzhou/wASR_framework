import click, numpy as np, ete3_extensions, os, subprocess, re, shutil, pandas as pd

treetime = '/home/zhemin/miniconda3/envs/p312/bin/treetime'
wASR = '/titan/softwares/EPHI/wASR/wASR.py'

def eval_treetime(confidence, ground_truth) :
    dat = pd.read_csv(confidence)
    res = [0, 0, 0, 0]
    for d in dat.values :
        if d[0] in ground_truth :
            trait, prob = ('x0', d[1]) if d[1] >= 0.5 else ('x1', d[2])
            i = 0 if prob >= 0.7 else 2
            i += 0 if trait == ground_truth[d[0]] else 1
            res[i] += 1
    return res


@click.command()
@click.option('-n', '--nexus')
def main(nexus) :
    results = {}
    # prepare folder
    odir = nexus.rsplit('.', 1)[0]
    if not os.path.isdir(odir) :
        os.makedirs(odir)
    
    # set up required files
    simtree = ete3_extensions.read_nexus(nexus)[0]
    simtree.name = 'NODE_0000000'
    ground_truth = {node.name:node.annotations['states'] for node in simtree.traverse() if not node.is_leaf()}
    traits = []
    for node in simtree.traverse() :
        node.dist *= 0.1
    
    for tip in simtree.get_leaves() :
        traits.append([tip.name, tip.annotations['states']])
    
    with open(os.path.join(odir, 'tree.nwk'), 'wt') as fout :
        fout.write(simtree.write(format=1))
    
    with open(os.path.join(odir, 'states.csv'), 'wt') as fout :
        fout.write('ID,state\n')
        for trait in traits :
            fout.write(f'{trait[0]},{trait[1]}\n')
    
    freq = [float(f) for f in re.findall(r'_x(\d+)_(\d+)_sim', nexus)[0]]

    # run wASR
    for err in np.arange(-30, 31, 2) :
        fx = 1+err/10. if err >= 0 else 1/(1-err/10)
        f2 = [freq[0]*fx, freq[1]]
        state_count = np.unique([t[1] for t in traits], return_counts=True)[1]
        base_num = 0 #np.log(len(traits))
        weights = [f2[0]/(state_count[0]+base_num), f2[1]/(state_count[1]+base_num)]

        with open(os.path.join(odir, 'states.weight'), 'wt') as fout :
            fout.write('state,weight\n')
            for i, w in enumerate(weights) :
                fout.write(f'x{i},{w}\n')

        subprocess.run(f'python {wASR} -w states.weight --states states.csv --tree tree.nwk --confidence --out iWeight_{err}.out'.split(), cwd=odir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        res = eval_treetime(f'{odir}/iWeight_{err}.out/confidence.csv', ground_truth)
        results[f'iWeight_{err}'] = res
 
    print(nexus, results)


if __name__ == '__main__' :
    main()
