import ete3_extensions, click, pandas as pd, numpy as np, os, collections


@click.command()
@click.option('-d', '--outdir', help='output folder')
@click.option('-n', '--nexus', help='nexus file')
@click.option('-t', '--trait', help='trait file')
@click.option('-w', '--weight', help='a number [default: 10] of maximum subsamples, or weight list of samples', default='10')
def main(outdir, nexus, trait, weight) :
    tre = ete3_extensions.read_nexus(nexus)[0]
    leaves = {n:'' for n in tre.get_leaf_names()}
    
    traits = collections.defaultdict(list)
    
    for n, t in pd.read_csv(trait, sep='\t', header=None, na_filter=False).values :
        if t and n in leaves :
            traits[t].append(n)
    
    subsamples = []
    try :
        weight = int(weight)
        for t, names in traits.items() :
            n = names if len(names) <= weight else np.random.choice(names, weight, replace=False).tolist()
            subsamples.extend([[nn, t] for nn in n])
    except :
        weights = pd.read_csv(weight, sep='\t', header=None, na_filter=False).values
        p = ['', 1e9]
        for t, w in weights[np.argsort(-weights.T[1])][:10] :
            n = len(traits.get(t, 0))
            if n > 1 and n/w < p[1] :
                p = [t, n/w]

        p_sample = {}
        for t, w in weights :
            n = len(traits.get(t, 0))
            p_sample[t] = max(w * p[1], 1.)

        for t, names in traits.items() :
            x = p_sample[t] # p_sample.get(t, 1)
            n_sample = int(x) + (1 if np.random.rand() < (x - int(x)) else 0)
            n = names if len(names) <= n_sample else np.random.choice(names, n_sample, replace=False).tolist()
            subsamples.extend([[nn, t] for nn in n])

    if not os.path.isdir(outdir) :
        os.makedirs(outdir)
    with open(os.path.join(outdir, 'dating.out.nex'), 'wt') as fout :
        fout.write(ete3_extensions.write_nexus([tre]))
    with open(os.path.join(outdir, 'trait.txt'), 'wt') as fout :
        fout.write('ID\tState\n')
        for n, t in subsamples :
            fout.write(f'{n}\t{t}\n')


if __name__ == '__main__' :
    main()
