import numpy as np, click, pandas as pd, gzip, collections, ete3_extensions, sys, ete3



def assign_mut_to_regions(mutation, dhmm, rhmm, node_conversion) :
    abbrs = { "Frameshift":"3_FS", "Nonsyn":"1_NS", 
                "Syn":"0_SY", "Indel":"4_ID",
                "Nonsense":"2_PM", "Intergenic":"6_IG", 
                "Upstream":"6_IG", "Undetermined":"5_UD" }

    diversifiedRegions = collections.defaultdict(list)
    dReg_tags = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(dhmm, 'rt') as fin :
        for line in fin :
            if line.startswith('\tDiversifiedRegion') :
                p = line.strip().split('\t')
                for i in range(int(p[3])//1000, 1+int(p[4])//1000) :
                    dReg_tags[p[2]][i].append(len(diversifiedRegions[p[2]]))
                diversifiedRegions[p[2]].append([int(p[3]), int(p[4])])
    
    importedRegions = collections.defaultdict(lambda: collections.defaultdict(list))
    rReg_tags = collections.defaultdict(lambda : collections.defaultdict(lambda: collections.defaultdict(list)))
    with open(rhmm, 'rt') as fin :
        for line in fin :
            if line.startswith('\tImportation') :
                p = line.strip().split('\t')
                for i in range(int(p[3])//1000, 1+int(p[4])//1000) :
                    rReg_tags[p[1]][p[2]][i].append(len(importedRegions[p[1]][p[2]]))
                importedRegions[p[1]][p[2]].append([int(p[3]), int(p[4])])
    
    branch_mutations = collections.defaultdict(lambda: collections.defaultdict(list))
    with gzip.open(mutation, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                continue
            p = line.strip().split('\t')
            branch, site = p[0], int(p[2])
            if branch not in node_conversion :
                continue
            branch = node_conversion[branch]
            
            muts = {}
            if len(p) > 5 :
                for g in p[5].split(';') :
                    gg = g.split(':')
                    if gg[2] != 'Upstream' :
                        muts[abbrs[gg[2]]] = 1
            else :
                muts['6_IG'] = 1
            
            gene_idx = 'normal'
            inRec = False
            if p[0] in importedRegions :
                for idx in rReg_tags[p[0]][p[1]][site//1000] :
                    s, e = importedRegions[p[0]][p[1]][idx]
                    if s <= site and site <= e :
                        inRec = True
                        break
            if not inRec :
                inDiv = False
                for idx in dReg_tags[p[1]][site//1000] :
                    s, e = diversifiedRegions[p[1]][idx]
                    if s <= site and site <= e :
                        inDiv = True
                        break
                if inDiv :
                    gene_idx = 'inDiv'
            else :
                gene_idx = 'inRec'
        
            mut_sum = sum(muts.values())
            for mut in muts.keys() :
                branch_mutations[(branch)][(gene_idx, mut)] += [(p[1], site)]
    return branch_mutations


@click.command()
@click.option('-p', '--prefix', help='prefix for outputs')
@click.option('-n', '--nexus', help='treetime nexus file')
@click.option('-R', '--root_node', help='root node for the tree')
@click.option('-t', '--tree', help='original tree for mutations')
@click.option('-m', '--mutations', help='mutations.gz file')
@click.option('-r', '--rhmm', help='rechmm file')
@click.option('-d', '--dhmm', help='divhmm file')
@click.option('-D', '--direction', default='down', type=click.Choice(['up', 'down']), help='direction of the tree')
def main(prefix, tree, mutations, rhmm, dhmm, nexus, direction, root_node) :
    tre = ete3.Tree(tree, format=1)
    if not tre.name :
        tre.name = 'NODE_0000000'

    timed_tre = ete3_extensions.read_nexus(nexus)[0]
    if not timed_tre.name :
        timed_tre.name = 'NODE_0000000'
    if root_node :
        timed_tre = [n for n in timed_tre.traverse('levelorder') if n.name == root_node][0]
        timed_tre.up = None
    
    scaling = 1
    timed_nodes = {timed_tre.name: [timed_tre, 0., 0.]}
    for n in timed_tre.get_descendants('preorder') :
        timed_nodes[n.name] = [n, timed_nodes[n.up.name][2], timed_nodes[n.up.name][2]+scaling*n.dist]

    paths = {}
    for tip in timed_tre.get_leaves() :
        paths[tip.name] = []
        for asc in [tip] + tip.get_ancestors() :
            if asc.up :
                if direction == 'up' :
                    paths[tip.name].append([asc.name, timed_nodes[asc.name][1], timed_nodes[asc.name][2]])
                else :
                    paths[tip.name].append([asc.name, timed_nodes[tip.name][2]-timed_nodes[asc.name][2], timed_nodes[tip.name][2]-timed_nodes[asc.name][1]])


    nodes = {}
    for n in tre.traverse() :
        if n.name in timed_nodes :
            nodes[n.name] = n.name
        else :
            for d in n.get_descendants('levelorder') :
                if d.name in timed_nodes :
                    nodes[n.name] = d.name
                    break

    branch_mutations = assign_mut_to_regions(mutations, dhmm, rhmm, nodes)
    
    timed_path = collections.defaultdict(set)
    for i, path in enumerate(paths.values()) :
        print('{0}/{1}          '.format(i+1, len(paths)), end='\r')
        for n, s, e in path :
            if e - s > 0.001 :
                for i in range(int(s+0.999999), int(e)+1, 1) :
                    # timed_path[i].append({t:{c:1/(e-s) for c in cs} for t, cs in branch_mutations.get(n, {}).items()})
                    # timed_path[i].append({t:len(cs)/(e-s) for t, cs in branch_mutations.get(n, {}).items()})
                    timed_path[i] |= set([n])
    
    res = {}
    for year, nodes in sorted(timed_path.items()) :
        print('{0}/{1}          '.format(year, len(timed_path)), end='\r')
        tp = [{t:len(cs)/(timed_nodes[n][2]-timed_nodes[n][1]) for t, cs in branch_mutations.get(n, {}).items()} for n in nodes] + [{}]
        dat = pd.DataFrame.from_dict(tp)
        dat.fillna(0, inplace=True)
        idx = np.random.choice(dat.shape[0], dat.shape[0]*3000)
        rdat = pd.DataFrame(np.sum(dat.values[idx].reshape([3000] + list(dat.shape)), 1)/(dat.shape[0]-1), columns=dat.columns)
        if ('normal', '0_SY') in rdat.columns :
            if ('inDiv', '0_SY') in rdat.columns :
                rdat[('inDiv-p', '0_SY')] = (rdat[('inDiv', '0_SY')]+1e-8)/(rdat[('normal', '0_SY')]+1e-8)
            if (('inRec-p', '0_SY')) in rdat.columns :
                rdat[('inRec-p', '0_SY')] = (rdat[('inRec', '0_SY')]+1e-8)/(rdat[('normal', '0_SY')]+1e-8)

        if ('normal', '1_NS') in rdat.columns :
            if ('inDiv', '1_NS') in rdat.columns :
                rdat[('inDiv-p', '1_NS')] = (rdat[('inDiv', '1_NS')]+1e-8)/(rdat[('normal', '1_NS')]+1e-8)
            if (('inRec-p', '1_NS')) in rdat.columns :
                rdat[('inRec-p', '1_NS')] = (rdat[('inRec', '1_NS')]+1e-8)/(rdat[('normal', '1_NS')]+1e-8)

        if ('normal', '2_PM') in rdat.columns :
            if ('inDiv', '2_PM') in rdat.columns :
                rdat[('inDiv-p', '2_PM')] = (rdat[('inDiv', '2_PM')]+1e-8)/(rdat[('normal', '2_PM')]+1e-8)
            if (('inRec-p', '2_PM')) in rdat.columns :
                rdat[('inRec-p', '2_PM')] = (rdat[('inRec', '2_PM')]+1e-8)/(rdat[('normal', '2_PM')]+1e-8)

        if ('normal', '3_FS') in rdat.columns :
            if ('inDiv', '3_FS') in rdat.columns :
                rdat[('inDiv-p', '3_FS')] = (rdat[('inDiv', '3_FS')]+1e-8)/(rdat[('normal', '3_FS')]+1e-8)
            if (('inRec-p', '3_FS')) in rdat.columns :
                rdat[('inRec-p', '3_FS')] = (rdat[('inRec', '3_FS')]+1e-8)/(rdat[('normal', '3_FS')]+1e-8)

        if ('normal', '4_ID') in rdat.columns :
            if ('inDiv', '4_ID') in rdat.columns :
                rdat[('inDiv-p', '4_ID')] = (rdat[('inDiv', '4_ID')]+1e-8)/(rdat[('normal', '4_ID')]+1e-8)
            if (('inRec-p', '4_ID')) in rdat.columns :
                rdat[('inRec-p', '4_ID')] = (rdat[('inRec', '4_ID')]+1e-8)/(rdat[('normal', '4_ID')]+1e-8)

        if ('normal', '6_IG') in rdat.columns :
            if ('inDiv', '6_IG') in rdat.columns :
                rdat[('inDiv-p', '6_IG')] = (rdat[('inDiv', '6_IG')]+1e-8)/(rdat[('normal', '6_IG')]+1e-8)
            if (('inRec-p', '6_IG')) in rdat.columns :
                rdat[('inRec-p', '6_IG')] = (rdat[('inRec', '6_IG')]+1e-8)/(rdat[('normal', '6_IG')]+1e-8)


        median, lower, upper = np.median(rdat, axis=0), np.quantile(rdat, 0.025, axis=0), np.quantile(rdat, 0.975, axis=0)
        
        res[year] = {'Num_path':len(tp)}
        for (r, t), m, l, u in zip(rdat.columns, median, lower, upper) :
            res[year][f'{r}|{t}|ave'] = m
            res[year][f'{r}|{t}|CI05'] = l
            res[year][f'{r}|{t}|CI95'] = u
        
    res = pd.DataFrame.from_dict(res)
    res = res.reindex(sorted(res.columns), axis=1).sort_index().T
    res2 = res.groupby(np.arange(len(res))//10).mean()
    res2.to_csv(prefix + '.timed_mutations.csv', sep=',')


if __name__ == '__main__' :
    main()
