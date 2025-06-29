#!/usr/bin/env python3
import numpy as np
import click
import pandas as pd
import gzip
import collections
import sys
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@click.command()
@click.option('-p', '--prefix', help='prefix for outputs', required=True)
@click.option('-g', '--geo_enrichment', help='gene_enrichment by DivHMM_to_trait', required=True)
@click.option('-k', '--gene_kos', help='gene_kos by DivHMM_to_genes', required=True)
def main(prefix, geo_enrichment, gene_kos):
    """
    Calculate KO enrichments for each enriched geographic region.
    """
    print(f"Loading geographic enrichment data from {geo_enrichment}")
    # Load the geographic enrichment data
    try:
        geo_df = pd.read_csv(geo_enrichment)
    except Exception as e:
        print(f"Error reading geo_enrichment file: {e}")
        sys.exit(1)
    
    # Filter regions where %pos > %neg and q-value < 0.05
    print("Filtering enriched regions...")
    enriched_regions = geo_df[(geo_df['%pos'] > geo_df['%neg']) & (geo_df['q_value'] < 0.05)]
    print(f"Found {len(enriched_regions)} enriched regions")
    
    # Load the gene_kos data
    print(f"Loading gene KO data from {gene_kos}")
    try:
        # Define column names based on the file format
        cols = ['gene_id', 'ko_id', 'contig', 'start', 'end', 'strand', 
                'snp_data1', 'snp_data2', 'snp_data3']
        ko_df = pd.read_csv(gene_kos, sep='\t', header=None, names=cols)
        
        diff_contig = ko_df['contig'][:-1].values != ko_df['contig'][1:].values
        
        prev_end = np.concatenate([[1], ko_df['end'][:-1]+1])
        prev_end[np.concatenate([[True], diff_contig])] = 1
        ko_df['pstart'] = [min(s1, s2) for s1, s2 in zip(ko_df['start'], prev_end)]
        post_start = np.concatenate([ko_df['start'][1:]-1, [9999999999]])
        post_start[np.concatenate([diff_contig, [True]])] = 9999999999
        ko_df['pend'] = [max(e1, e2) for e1, e2 in zip(ko_df['end'], post_start)]
        
    except Exception as e:
        print(f"Error reading gene_kos file: {e}")
        sys.exit(1)
    
    # Extract regions and their coordinates from enriched_regions
    region_coords = {}
    for _, row in enriched_regions.iterrows():
        # Check if it's an inDiv region (contains coordinates)
        if 'inDiv:' in str(row['Region']):
            if row['Region'] in region_coords :
                region_coords[row['Region']]['geo_region'].append(row['Trait_group'])
            else :
                # Parse region string to extract coordinates
                region_parts = row['Region'].split(':')
                if len(region_parts) >= 2:
                    contig_coord = region_parts[1].rsplit('_', 2)
                    if len(contig_coord) >= 3:
                        contig = contig_coord[0]
                        start = int(contig_coord[1])
                        end = int(contig_coord[2])
                        region_coords[row['Region']] = {'contig': contig, 'start': start, 'end': end, 
                                                        'region_name': row['Region'],
                                                        'geo_region': [row['Trait_group']]}
    
    print(f"Extracted coordinates for {len(region_coords)} enriched regions")
    
    # Find KOs that overlap with each region
    region_ko_map = collections.defaultdict(list)
    for _, gene in ko_df.iterrows():
        gene_contig = gene['contig']
        gene_start = int(gene['pstart'])
        gene_end = int(gene['pend'])
        
        # Check overlap with all regions
        for region_id, region in region_coords.items():
            if (gene_contig == region['contig'] and 
                max(gene_start, region['start']) <= min(gene_end, region['end'])):
                region_ko_map[region_id].append((gene['gene_id'], gene['ko_id']))
    
    region_genes = collections.defaultdict(set)
    
    for region_id, region_info in region_coords.items() :
        for geo_region in region_info['geo_region'] :
            region_genes[geo_region] |= set(region_ko_map[region_id])
    
    with open(f'{prefix}.DHMM.geo_gene.enrichment', 'wt') as fout :
        for geo_region, genes in sorted(region_genes.items()) :
            for gene, ko in sorted(genes) :
                fout.write(f'{gene}\t{ko}\t{geo_region}\n')
    
    

if __name__ == '__main__':
    main()