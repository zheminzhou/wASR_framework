import numpy as np, click, pandas as pd, gzip, collections, ete3_extensions, sys, ete3
from scipy.stats import fisher_exact, false_discovery_control


def assign_mut_to_regions(mutation, dhmm, rhmm, node_conversion) :
    abbrs = { "Frameshift":"3_FS", "Nonsyn":"1_NS", 
                "Syn":"0_SY",      "Indel":"4_ID",
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
                inDiv = None
                for idx in dReg_tags[p[1]][site//1000] :
                    s, e = diversifiedRegions[p[1]][idx]
                    if s <= site and site <= e :
                        inDiv = f'inDiv:{p[1]}_{s}_{e}'
                        break
                if inDiv :
                    gene_idx = inDiv
            else :
                gene_idx = 'inRec'
        
            # mut_sum = sum(muts.values())
            for mut in muts.keys() :
                branch_mutations[(branch)][gene_idx] += [mut]
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
    global trait_conv
    tre = ete3.Tree(tree, format=1)
    if not tre.name :
        tre.name = 'NODE_0000000'

    timed_tre = ete3_extensions.read_nexus(nexus)[0]
    if not timed_tre.name :
        timed_tre.name = 'NODE_0000000'
    if root_node :
        timed_tre = [n for n in timed_tre.traverse('levelorder') if n.name == root_node][0]
        timed_tre.up = None
        
    timed_nodes = {}
    for n in timed_tre.traverse('preorder') :
        timed_nodes[n.name] = [n, trait_conv.get(n.annotations['state'], n.annotations['state']), n.annotations['state.prop']]

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

    trait_conv = {}
    trait_mutations = collections.defaultdict(lambda: collections.defaultdict(dict))
    for name, (node, trait, prop) in timed_nodes.items() :
        for reg, mut in branch_mutations.get(name, {}).items() :
            trait_mutations[reg][trait_conv.get(trait, trait)][name] = len(mut)
    
    res = {}
    
    for region, mutations in trait_mutations.items() :
        dat = pd.DataFrame.from_dict(mutations)
        dat.fillna(0., inplace=True)
        # idx = np.random.choice(dat.shape[0], dat.shape[0]*3000)
        # rdat = pd.DataFrame(np.sum(dat.values[idx].reshape([3000] + list(dat.shape)), 1), columns=dat.columns)
        res[region] = dat.sum()
    res = pd.DataFrame(res).fillna(0.0).T
    for r, n in res.sum(axis=0).items() :
        print(f'{r}\t{n}')
    res.to_csv(prefix + '.trait_mutations.csv', sep=',')
    
    tests = []
    negatives = res.loc['normal']
    for reg, positives in res.iterrows() :
        if reg == 'normal' : continue
        n_neg = negatives.sum()
        n_pos = positives.sum()
        # if n_pos < 3 :
        #     continue
        data = list(zip(positives.index, positives, negatives))
        # data = [d for d in data if d[0] != 'ND']

        # for t, p, n in data :
        #     pvalues = fisher_exact([[p, n_pos - p], [n, n_neg - n]])
        #     tests.append([reg, t, p, p/(n_pos + 1e-8), n/(n_neg + 1e-8), pvalues[1]])

        while len(data) :
            test = []
            pprev = 0
            for t, p, n in data :
                pvalues = fisher_exact([[p, n_pos - p], [n, n_neg - n]])
                test.append([reg, t, p, p/(n_pos + 1e-8), n/(n_neg + 1e-8), pvalues[1]])
            
            idx = np.argmin([t[-1] for t in test])
            if test[idx][5] < pprev :
                test[idx][5] = pprev
            else :
                pprev = test[idx][5]
            tests.append(test[idx])
            n_neg, n_pos = n_neg - data[idx][2], n_pos - data[idx][1]
            data = data[:idx] + data[idx+1:]

    fdrs = false_discovery_control([t[-1] for t in tests])
    for fdr, t in zip(fdrs, tests) :
        t.append(fdr)
    tests = pd.DataFrame(tests, columns=['Region', 'Trait_group', 'n_SNP', '%pos', '%neg', 'p_value', 'q_value'])
    tests = tests.sort_values(by=['p_value'])
    tests = tests.loc[tests['p_value'] < 0.05]
    tests.to_csv(prefix + '.trait_mutations.enrichment', sep=',')


Bp_conv = '''United_States	United_States
France	France
Japan	Japan
China	China
Australia	Oceania
United_Kingdom	United_Kingdom
Kenya	Kenya
Netherlands	Netherlands
Spain	Spain
Norway	Norway
Austria	Austria
Tunisia	Africa_North
Brazil	Americas
Argentina	Americas
Canada	Americas
Mexico	Americas
Guatemala	Americas
Haiti	Americas
Iran	Asia
India	Asia
Israel	Asia
Vietnam	Asia
Czech_Republic	Europe
Sweden	Europe
Italy	Europe
Ireland	Europe
Finland	Europe
Russia	Europe
Poland	Europe
Belgium	Europe
Denmark	Europe
New_Zealand	Oceania'''
Bp_conv = { line.split('\t')[0]:line.split('\t')[1] for line in Bp_conv.split('\n') }


ParaA_conv = '''Morocco	Africa
Senegal	Africa
South_Africa	Africa
Sudan	Africa
Tunisia	Africa
Mali	Africa
Algeria	Africa
Malawi	Africa
Sierra_Leone	Africa
Guinea	Africa
Gambia	Africa
Ghana	Africa
Ethiopia	Africa
Chad	Africa
Zambia	Africa
Brazil	Americas
Peru	Americas
Argentina	Americas
Australia	Australia
Japan	Asia_East
Singapore	Asia_SE
Cambodia	Asia_SE
Indonesia	Asia_SE
Vietnam	Asia_SE
Myanmar	Asia_SE
Sri_Lanka	Sri_Lanka
Bangladesh	Bangladesh
China_guandong	China
China_zhejiang	China
China_jiangsu	China
China_jiangxi	China
China	China
China_beijing	China_N
China_henan	China_N
China_shandong	China_N
China_shanghai	China
China_hunan	China
Taiwan	China
China_sichuan	China
China_fujian	China
China_xinjiang	China_N
China_tianjin	China_N
China_hongkong	China
China_guangxi	China_SW
China_yunnan	China_SW
China_guizhou	China_SW
France	Europe
Portugal	Europe
Ireland	Europe
Norway	Europe
Albania	Europe
Bulgaria	Europe
Denmark	Europe
Luxembourg	Europe
India	India
ND	ND
Turkey	Asia_West
Jordan	Asia_West
Kuwait	Asia_West
Malta	Europe
Lebanon	Asia_West
Palestine	Asia_West
Iran	Asia_West
Nepal	Nepal
United_States	North_Americas
Canada	North_Americas
Pakistan	Pakistan
United_Kingdom	United_Kingdom'''
ParaA_conv = { line.split('\t')[0]:line.split('\t')[1] for line in ParaA_conv.split('\n') }


TB_conv = '''Malawi	Africa_East
Madagascar	Africa_East
Uganda	Africa_East
Tanzania	Africa_East
DR_Congo	Africa_East
Djibouti	Africa_East
Mozambique	Africa_East
Kenya	Africa_East
Zimbabwe	Africa_East
Ethiopia	Africa_East
Tunisia	Africa_North
Morocco	Africa_North
South_Africa	Africa_South
Botswana	Africa_South
Namibia	Africa_South
Gambia	Africa_West
Nigeria	Africa_West
Ghana	Africa_West
Liberia	Africa_West
Cote_d'Ivoire	Africa_West
Mali	Africa_West
Sierra_Leone	Africa_West
Kazakhstan	Asia_Central
Tajikistan	Asia_Central
Kyrgyzstan	Asia_Central
Japan	Asia_East
South_Korea	Asia_East
Thailand	Asia_SE
Indonesia	Asia_SE
Malaysia	Asia_SE
Vietnam	Asia_SE
Cambodia	Asia_SE
Myanmar	Asia_SE
Philippines	Asia_SE
Bangladesh	Asia_South
Pakistan	Asia_South
Nepal	Asia_South
Sri_Lanka	Asia_South
Bhutan	Asia_South
Georgia	Asia_West
Iran	Asia_West
Azerbaijan	Asia_West
Israel	Asia_West
Lebanon	Asia_West
Oman	Asia_West
Turkey	Asia_West
China	Asia_East
Russia	Europe_East
Moldova	Europe_East
Romania	Europe_East
Belarus	Europe_East
Ukraine	Europe_East
Czechia	Europe_East
Slovakia	Europe_East
Bulgaria	Europe_East
Hungary	Europe_East
Poland	Europe_East
Sweden	Europe_North
Ireland	Europe_North
Denmark	Europe_North
Norway	Europe_North
Finland	Europe_North
Estonia	Europe_North
Italy	Europe_South
Portugal	Europe_South
Spain	Europe_South
Albania	Europe_South
Serbia	Europe_South
Bosnia_and_Herzegovina	Europe_South
Croatia	Europe_South
Switzerland	Europe_West
Netherlands	Europe_West
France	Europe_West
Belgium	Europe_West
Germany	Europe_West
India	Asia_South
United_States	North_America
Canada	North_America
Australia	Oceania
Papua_New_Guinea	Oceania
New_Zealand	Oceania
Brazil	South_America
Peru	South_America
Mexico	South_America
Colombia	South_America
Panama	South_America
Guatemala	South_America
Paraguay	South_America
Argentina	South_America
Uruguay	South_America
Ecuador	South_America
United_Kingdom	Europe_West'''
TB_conv = { line.split('\t')[0]:line.split('\t')[1] for line in TB_conv.split('\n') }

trait_conv = {}
trait_conv = TB_conv


if __name__ == '__main__' :
    main()
