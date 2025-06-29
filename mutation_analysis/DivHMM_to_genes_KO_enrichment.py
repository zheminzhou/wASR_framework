import os, click, collections, gzip, re, json, numpy as np, pandas as pd
from scipy.stats import fisher_exact, false_discovery_control
from keggtools import Enrichment, Resolver



brite_groups = []
def iter_brites(brites) :
    if 'children' in brites :
        tmp = []
        for child in brites['children'] :
            tmp.extend(iter_brites(child))
        brite_groups.append([brites['name'], tmp])
        return tmp
    else :
        return [brites['name'].split()[0]]
        

@click.command()
@click.option('-o', '--outfile', help='outputs')
@click.option('-d', '--dhmm_genes', help='*.diversified.related_genes')
@click.option('-k', '--ko_brite', help='brite file', default='/titan/softwares/RecHMM/ko00001.json', show_default=True)
def ko_enrichment(outfile, dhmm_genes, ko_brite) :
    if not outfile :
        outfile = os.path.basename(dhmm_genes) + '.enrichment'

    snp_categories = [[], [], []]
    gene_categories = [[[], []], [[], []]]
    with open(dhmm_genes, 'rt') as fin :
        for line in fin :
            p = line.strip().split('\t')
            snp_cats = [json.loads(json_str) for json_str in p[6:9]]
            for sc, x in zip(snp_categories, snp_cats) :
                sc.append(x)

            if p[1] != '-' :
                x = {k:v for k, v in snp_cats[1].items()}
                if len(x) > 0 :
                    gene_categories[0][1].append(p[1])
                else :
                    gene_categories[0][0].append(p[1])
                
                if len(snp_cats[2]) > 0 :
                    gene_categories[1][1].append(p[1])
                else :
                    gene_categories[1][0].append(p[1])

    snp_res = [{}, {}, {}]

    for idx in range(3) :
        dat = pd.DataFrame.from_dict(snp_categories[idx])
        dat.fillna(0, inplace=True)
        ridx = np.random.choice(dat.shape[0], dat.shape[0]*1000)
        rdat = pd.DataFrame(np.sum(dat.values[ridx].reshape([1000, dat.shape[0], dat.shape[1]]), (1)), columns=dat.columns)
        rdat['7_dNdS'] = (rdat['1_NS'] + rdat['2_PM'])/rdat['0_SY']
        rdat['8_dPseu'] = (rdat['2_PM'] + rdat['3_FS'])/rdat['0_SY']
        rdat['9_dIntergenic'] = rdat['6_IG']/rdat['0_SY']
        lower = np.quantile(rdat, 0.025, axis=0)
        upper = np.quantile(rdat, 0.975, axis=0)
        mean = np.mean(rdat, 0)
        for c, m, l, u in zip(rdat.columns, mean, lower, upper) :
            snp_res[idx][c] = f'{m:.3f}, {l:.3f}, {u:.3f}'
    
    with open(outfile, 'wt') as fout :
        print(json.dumps(snp_res, indent=2, sort_keys=True))
        fout.write(json.dumps(snp_res, sort_keys=True, indent=2)+'\n')
        
        brites = json.load(open(ko_brite, 'rt'))
        iter_brites(brites)

        brite_info = []
        for tag, brite in brite_groups :
            groups = [[len(set(brite) & set(gene_categories[0][0])), len(set(brite) & set(gene_categories[0][1]))],
                    [len(set(brite) & set(gene_categories[1][0])), len(set(brite) & set(gene_categories[1][1]))]]
            brite_info.append([tag] + groups)
        
        base_line = brite_info[-1][1]
        results = []
        for tag, groups, g2 in brite_info[:-1] :
            if sum(groups) <= 3 or groups[0] <= 2 :
                continue
            tests = [[groups[0], groups[1]], [base_line[0]-groups[0], base_line[1]-groups[1]]]
            res = fisher_exact(tests, alternative='less')
            p0 = groups[1]/(groups[0]+groups[1])
            p1 = base_line[1]/(base_line[0]+base_line[1])
            results.append([tag.replace(' ', '_'), groups[1], groups[0], p0, p1, res[1]])
        
        paths  = [r for r in results if r[0].find('PATH') >= 0]
        brites = [r for r in results if r[0].find('BR:') >= 0]
        groups = [r for r in results if r[0].find(':ko') < 0]
        
        fdr1 = false_discovery_control([r[5] for r in paths])
        fdr2 = false_discovery_control([r[5] for r in brites])
        fdr3 = false_discovery_control([r[5] for r in groups])
        
        results = paths + brites + groups
        fdr = np.concatenate([fdr1 , fdr2 , fdr3])
        
        for idx in np.argsort(fdr) :
            if results[idx][-1] < 0.05 :
                print(f'Div - {results[idx][0]}:\t{results[idx][1]}\t{results[idx][2]}\t{results[idx][3]:.3f}\t{results[idx][4]:.3f}\t{results[idx][5]:.3e}\t{fdr[idx]:.3e}')
                fout.write(f'Div - {results[idx][0]}:\t{results[idx][1]}\t{results[idx][2]}\t{results[idx][3]:.3f}\t{results[idx][4]:.3f}\t{results[idx][5]:.3e}\t{fdr[idx]:.3e}\n')


        base_line = brite_info[-1][2]
        results = []
        for tag, g1, groups in brite_info[:-1] :
            if sum(groups) <= 3 :
                continue
            tests = [[groups[0], groups[1]], [base_line[0]-groups[0], base_line[1]-groups[1]]]
            res = fisher_exact(tests, alternative='less')
            p0 = groups[1]/(groups[0]+groups[1])
            p1 = base_line[1]/(base_line[0]+base_line[1])
            results.append([tag.replace(' ', '_'), groups[1], groups[0], p0, p1, res[1]])
        
        paths  = [r for r in results if r[0].find('PATH') >= 0]
        brites = [r for r in results if r[0].find('BR:') >= 0]
        groups = [r for r in results if r[0].find(':ko') < 0]
        
        fdr1 = false_discovery_control([r[5] for r in paths])
        fdr2 = false_discovery_control([r[5] for r in brites])
        fdr3 = false_discovery_control([r[5] for r in groups])
        
        results = paths + brites + groups
        fdr = np.concatenate([fdr1 , fdr2 , fdr3])
        
        for idx in np.argsort(fdr) :
            if results[idx][-1] < 0.05 :
                print(f'Rec - {results[idx][0]}:\t{results[idx][1]}\t{results[idx][2]}\t{results[idx][3]:.3f}\t{results[idx][4]:.3f}\t{results[idx][5]:.3e}\t{fdr[idx]:.3e}')
                fout.write(f'Rec - {results[idx][0]}:\t{results[idx][1]}\t{results[idx][2]}\t{results[idx][3]:.3f}\t{results[idx][4]:.3f}\t{results[idx][5]:.3e}\t{fdr[idx]:.3e}\n')


if __name__ == '__main__' :
    ko_enrichment()
