import ete3_extensions, click, ete3, numpy as np


@click.command()
@click.argument('nexus')
def main(nexus) :
    prefix = nexus.rsplit('.', 1)[0]
    with open(nexus, 'rt') as fin :
        for line in fin :
            if line.startswith('BEGIN TREES;') :
                tre_line = fin.readline().split('=', 1)[-1].split(';')[0].strip() +';';
                tre = ete3.Tree(tre_line, format=1)
            elif line.startswith('  MATRIX') :
                break
        node_info = {}
        for line in fin :
            if line.startswith('  ;') :
                break
            p = line.strip().split()
            node_info[p[0]] = int(p[1])-1
    states = [[], []]
    for n in tre.traverse('postorder') :
        s = node_info[n.name]
        n.annotations = {'states':f'x{s}'}
        if n.is_leaf() :
            states[s].append(n.name)

    patterns = [[1.,1], [1.,3], [1.,10], [1.,30]]
    for p in patterns :
        for ite in np.arange(3) :
            n1 = min(len(states[1]), 200)
            n0 = n1 * p[0] / p[1]
            n0 = max(min(int(n0) + (np.random.rand() < n0-int(n0)), len(states[0])), 1)
            rtip = np.random.choice(states[0], n0, replace=False).tolist() + np.random.choice(states[1], n1, replace=False).tolist()
            t0 = tre.copy()
            t0 = ete3_extensions.prune(t0, rtip)
            outfile = f'{prefix}.sub_{int(p[0])}_{p[1]}.ite_{ite}.nex'
            with open(outfile, 'wt') as fout :
                fout.write(ete3_extensions.write_nexus([t0]))

if __name__ == '__main__' :
    main()
