import click
import pandas as pd
import math
from collections import defaultdict

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file path. If not provided, prints to stdout')
def calculate_similarity(input_file, output):
    """Calculate the proportion of shared blocks between different regions using JC distance."""
    # Read the input CSV file
    df = pd.read_csv(input_file)
    
    # Extract regions and their blocks
    region_blocks = defaultdict(set)
    all_blocks = set()
    
    for _, row in df.iterrows():
        if row['q_value'] >= 0.01 :
            continue
        if row['%pos'] < row['%neg'] :
            continue
        region = row['Trait_group']
        block = row['Region']  # Second column contains block information
        region_blocks[region].add(block)
        all_blocks.add(block)

    # Find all unique regions
    regions = list(region_blocks.keys())
    
    # Calculate Jukes-Cantor distance and similarity for each pair of regions
    results = []
    
    for i in range(len(regions)):
        for j in range(i+1, len(regions)):
            region1 = regions[i]
            region2 = regions[j]
            
            blocks1 = region_blocks[region1]
            blocks2 = region_blocks[region2]
            
            # Calculate proportion of shared blocks
            shared_blocks = blocks1.intersection(blocks2)
            total_unique_blocks = blocks1.union(blocks2)
            similarity = len(shared_blocks)/(len(total_unique_blocks)+5)
            
            results.append((region1, region2, similarity))
    
    # Sort results by similarity (descending)
    results.sort(key=lambda x: x[2], reverse=True)
    
    # Prepare the output
    output_lines = ["country1,country2,similarity"]
    for region1, region2, similarity in results:
        if similarity > 0 :
            output_lines.append(f"{region1},{region2},{similarity:.6f}")
    
    output_text = "\n".join(output_lines)
    
    # Write output to file or print to stdout
    if output:
        with open(output, 'w') as f:
            f.write(output_text)
        click.echo(f"Results written to {output}")
    else:
        click.echo(output_text)

if __name__ == '__main__':
    calculate_similarity()
