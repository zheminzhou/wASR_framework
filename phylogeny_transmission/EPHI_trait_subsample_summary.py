import ete3_extensions, click, pandas as pd, numpy as np, os, collections


@click.command()
@click.option('-o', '--out_nexus', help='nexus file')
@click.option('-l', '--list', help='list of subsamples')
def main(out_nexus, list) :
    node_status = collections.defaultdict(lambda : collections.defaultdict(float))
    with open(list, 'rt') as fin :
        for line in fin :
            fn = line.strip().split()[0]
            tre = ete3_extensions.read_nexus(fn)[0]
            for n in tre.traverse() :
                idx = np.argmax(n.annotations['traits.prop'])
                trait, prop = n.annotations['traits.list'][idx], n.annotations['traits.prop'][idx]
                node_status[n.name][trait] += prop
    for n in tre.traverse() :
        traits, props = zip(*node_status[n.name].items())
        idx = np.argmax(props)
        max_state, max_prop = traits[idx], props[idx]/sum(props)
        n.annotations['state'], n.annotations['state.prop'] = max_state, max_prop
        n.annotations['traits.list'], n.annotations['traits.prop'] = traits, props
    
    with open(out_nexus, 'wt') as fout :
        fout.write(ete3_extensions.write_nexus([tre]))


if __name__ == '__main__' :
    main()
