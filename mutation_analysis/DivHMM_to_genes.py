import os, click, collections, gzip, re, json


@click.command()
@click.option('-o', '--outfile', help='outputs')
@click.option('-d', '--dhmm', help='*.diversified.region')
@click.option('-r', '--rhmm', help='*.recombination.region')
@click.option('-m', '--mutation', help='*.mutations.gz')
@click.option('-g', '--gff', help='gff file for gene annotations')
@click.option('-t', '--gene_tag', help='tag for genes', default='locus_tag', show_default=True)
@click.option('-k', '--ko_list', help='ko file for genes')
def get_genes(outfile, dhmm, rhmm, mutation, gff, gene_tag, ko_list) :
    if not outfile :
        outfile = os.path.basename(dhmm) + '.related_genes'

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
    
    genes = {}
    genomes = collections.defaultdict(lambda: collections.defaultdict(list))
    with gzip.open(gff, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                continue
            p = line.strip().split('\t')
            if p[2] in ('gene', 'pseudogene'):
                gene_name = re.findall(r'{0}=([^;]+)'.format(gene_tag), p[8])[0]
                s, e = int(p[3]), int(p[4])
                genes[gene_name] = [p[0], '-', s, e, p[6], {}, {}, {}]
                for i in range(s, e+1) :
                    genomes[p[0]][i].append(gene_name)
    
    with open(ko_list, 'rt') as fin :
        for line in fin :
            p = line.strip().split()
            if len(p) > 1 :
                gene, ko = p[0].split(':', 1)[-1], p[1].split(':', 1)[-1]
                if gene in genes :
                    genes[gene][1] = ko
    
    with gzip.open(mutation, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                if line.startswith('## Sequence_length:') :
                    p = line.strip().split()
                    seqLen = int(p[3])
                    prev = []
                    for i in range(1, seqLen+1) :
                        if genomes[p[2]][i] :
                            prev = genomes[p[2]][i]
                        else :
                            genomes[p[2]][i] = ['I'] + prev
                    prev = []
                    for i in range(seqLen, 0, -1) :
                        if genomes[p[2]][i][0] != 'I' :
                            prev = genomes[p[2]][i]
                        else :
                            genomes[p[2]][i].extend(prev)
                    
                # elif line.startswith('## Missing_region:') :
                #     p = line.strip().split()
                #     for i in range(int(p[3]), int(p[4])+1) :
                #         if genomes[p[2]][i][0] == 'I' :
                #             genomes[p[2]][i][0] = '-'
                #         else :
                #             genomes[p[2]][i] = ['-'] + genomes[p[2]][i]
                continue
            p = line.strip().split('\t')
            site = int(p[2])
            regions = genomes[p[1]][site]
            abbrs = {"Frameshift":"3_FS", "Nonsyn":"1_NS", 
                     "Syn":"0_SY", "Indel":"4_ID", 
                     "Nonsense":"2_PM", "Intergenic":"6_IG", 
                     "Undetermined":"5_UD"}
            muts = {r:'6_IG' for r in regions[1:]} if regions[0] == 'I' else {r:'5_UD' for r in regions[1:]}
            if len(p) > 5 :
                for g in p[5].split(';') :
                    gg = g.split(':')
                    if gg[2] != 'Upstream' :
                        muts[gg[0]] = abbrs[gg[2]]
            
            gene_idx = 5
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
                    gene_idx = 6
            else :
                gene_idx = 7
        
            for gene, mut in muts.items() :
                genes[gene][gene_idx][mut] = genes[gene][gene_idx].get(mut, 0) + 1
    
    with open(outfile, 'wt') as fout :
        for gene, (cont, ko, s, e, d, t1, t2, t3) in sorted(genes.items()) :
            t1, t2, t3 = json.dumps(t1, sort_keys=True).replace(' ', ''), \
                         json.dumps(t2, sort_keys=True).replace(' ', ''), \
                         json.dumps(t3, sort_keys=True).replace(' ', '')
            fout.write(f'{gene}\t{ko}\t{cont}\t{s}\t{e}\t{d}\t{t1}\t{t2}\t{t3}\n')



if __name__ == '__main__' :
    get_genes()
