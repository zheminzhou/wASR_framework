import os, click, collections, gzip, re, json


@click.command()
@click.option('-o', '--outfile', help='outputs')
@click.option('-d', '--dhmm', help='*.diversified.region')
@click.option('-r', '--rhmm', help='*.recombination.region')
@click.option('-m', '--mutation', help='*.mutations.gz')
@click.option('-g', '--gff', help='gff file for gene annotations')
@click.option('-t', '--gene_tag', help='tag for genes', default='locus_tag', show_default=True)
def get_genes(outfile, dhmm, rhmm, mutation, gff, gene_tag) :
    if not outfile :
        outfile = os.path.basename(dhmm) + '.annotations'

    diversifiedRegions = collections.defaultdict(list)
    dReg_tags = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(dhmm, 'rt') as fin :
        for line in fin :
            if line.startswith('\tDiversifiedRegion') :
                p = line.strip().split('\t')
                for i in range(int(p[3])//1000, 1+int(p[4])//1000) :
                    dReg_tags[p[2]][i].append(len(diversifiedRegions[p[2]]))
                diversifiedRegions[p[2]].append([int(p[3]), int(p[4]), [0, 0, set([])]])
    
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
                    
                continue
            p = line.strip().split('\t')
            site = int(p[2])
            regions = genomes[p[1]][site]
            mtype = 0 if len(p[4]) <= 4 else 1
            if regions[0] == 'I' :
                associated_genes = {f'I_{r}' for r in regions[1:]}
            else :
                associated_genes = set(regions)
            
            inRec = False
            if p[0] in importedRegions :
                for idx in rReg_tags[p[0]][p[1]][site//1000] :
                    s, e = importedRegions[p[0]][p[1]][idx]
                    if s <= site and site <= e :
                        inRec = True
                        break
            if not inRec :
                for idx in dReg_tags[p[1]][site//1000] :
                    s, e, m = diversifiedRegions[p[1]][idx]
                    if s <= site and site <= e :
                        m[mtype] += 1
                        m[2] |= associated_genes
    
    with open(outfile, 'wt') as fout :
        for cont, region in diversifiedRegions.items() :
            for s, e, (mut, indel, genes) in region :
                if all([g.startswith('I') for g in genes]) :
                    genes = ",".join(sorted(genes))
                else :
                    genes = ",".join(sorted([g for g in genes if not g.startswith('I_')]))
                fout.write(f'{cont}\t{s}\t{e}\t{mut}\t{indel}\t{genes}\n')

    print(f"Genes and mutations written to {outfile}")


if __name__ == '__main__' :
    get_genes()
