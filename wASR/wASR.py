import os, sys
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio import __version__ as bioversion
from treeanc import TreeAnc
from gtr import GTR
import click
from datetime import datetime
import collections


def get_outdir(outdir, suffix='_wASR'):
    if outdir:
        if os.path.exists(outdir):
            return outdir.rstrip('/') + '/' if os.path.isdir(outdir) else None
        else:
            os.makedirs(outdir)
            return outdir.rstrip('/') + '/'
    
    outdir_stem = datetime.now().date().isoformat()
    outdir = outdir_stem + suffix.rstrip('/') + '/'
    count = 1
    while os.path.exists(outdir):
        outdir = f"{outdir_stem}-{count:04d}{suffix.rstrip('/')}/"
        count += 1
    
    os.makedirs(outdir)
    return outdir


def parse_weights(leaf_to_attributes, frequencies=None, weights=None, power=0.5):
    if not weights:
        if frequencies:
            attribute_counts = collections.Counter(leaf_to_attributes.values())
            base_count = np.log(len(leaf_to_attributes.keys()))
            try:
                tmp_frequencies = pd.read_csv(frequencies, sep='\t' if frequencies.endswith('tsv') else ',', skipinitialspace=True)
                weight_dict = {state: float(freq) / (base_count + attribute_counts.get(state, 0)) 
                              for state, freq in tmp_frequencies.values if isinstance(freq, (int, float, str))}
            except:
                raise ValueError(f"Loading of frequencies file '{frequencies}' failed!")
        else:
            weight_dict = {state: 1.0 for state in leaf_to_attributes.values()}
    else:
        try:
            tmp_weights = pd.read_csv(weights, sep='\t' if weights.endswith('tsv') else ',', skipinitialspace=True)
            weight_dict = {state: float(weight) for state, weight in tmp_weights.values 
                          if isinstance(weight, (int, float, str))}
        except:
            raise ValueError(f"Loading of weights file '{weights}' failed!")
    
    weight_dict = {state: weight**power for state, weight in weight_dict.items()}
    mean_weight = np.mean([weight_dict[state] for state in leaf_to_attributes.values() if state in weight_dict])
    return {leaf: weight_dict.get(state, mean_weight) / mean_weight for leaf, state in leaf_to_attributes.items()}


def reconstruct_discrete_traits(tree, traits, missing_data='?', pc=1.0, sampling_bias_correction=None,
                                weights=None, verbose=0, iterations=7, rng_seed=None):
    
    unique_states = sorted(set(traits.values()) - {missing_data})
    n_states = len(unique_states)
    
    if n_states < 2:
        print("mugration: only one or zero states found -- this doesn't make any sense", file=sys.stderr)
        return None, None, None
    
    # Create alphabet mapping
    reverse_alphabet = {state: chr(65 + i) for i, state in enumerate(unique_states)}
    alphabet = list(reverse_alphabet.values())
    letter_to_state = {v: k for k, v in reverse_alphabet.items()}
    
    missing_char = chr(65 + n_states)
    reverse_alphabet[missing_data] = missing_char
    letter_to_state[missing_char] = missing_data
    
    # Construct GTR model
    W = np.ones((n_states, n_states), dtype=float)
    mugration_GTR = GTR.custom(pi=None, W=W, alphabet=np.array(alphabet))
    mugration_GTR.profile_map[missing_char] = np.ones(n_states)
    mugration_GTR.ambiguous = missing_char
    
    # Set up TreeAnc
    treeanc = TreeAnc(tree, gtr=mugration_GTR, verbose=verbose, ref='A',
                      convert_upper=False, one_mutation=0.0001, rng_seed=rng_seed)
    treeanc.use_mutation_length = False
    
    pseudo_seqs = {n.name: {0: reverse_alphabet[traits[n.name]] if n.name in traits else missing_char}
                   for n in treeanc.tree.get_terminals()}
    
    valid_seq = np.array([s[0] != missing_char for s in pseudo_seqs.values()])
    print(f"Assigned discrete traits to {np.sum(valid_seq)} out of {len(valid_seq)} taxa.\n")
    
    treeanc.aln = pseudo_seqs
    for node in treeanc.tree.get_terminals():
        node.weight = weights[node.name] if weights else 1.0
    
    try:
        treeanc.infer_ancestral_sequences(method='ml', infer_gtr=True, store_compressed=False, 
                                         pc=pc, marginal=True, normalized_rate=False, 
                                         fixed_pi=None, reconstruct_tip_states=True)
        treeanc.optimize_gtr_rate()
    except:
        raise Exception("Ancestral reconstruction failed, please see above for error messages and/or rerun with --verbose 4")
    
    for i in range(iterations):
        treeanc.infer_gtr(marginal=True, normalized_rate=False, pc=pc, fixed_pi=None)
        treeanc.optimize_gtr_rate()
    
    if sampling_bias_correction:
        treeanc.gtr.mu *= sampling_bias_correction
    
    treeanc.infer_ancestral_sequences(infer_gtr=False, store_compressed=False,
                                     marginal=True, normalized_rate=False,
                                     reconstruct_tip_states=True)
    
    return treeanc, letter_to_state, reverse_alphabet


@click.command()
@click.option('-t', '--tree', required=True, type=str, help="Tree file in newick, nexus, or phylip format")
@click.option('-s', '--states', required=True, type=str, help="CSV/TSV file with discrete characters")
@click.option('-f', '--frequencies', help="CSV/TSV file with expected total counts of each state")
@click.option('-w', '--weights', help="CSV/TSV file with expected sampling rate of each state")
@click.option('-p', '--power', help="Power adjustment of the weights", default=1.0)
@click.option('-r', '--rng_seed', type=int, help="Random number generator seed for treetime")
@click.option('-n', '--name_column', help="Label of the column to be used as taxon name")
@click.option('-a', '--attribute', help="Attribute to reconstruct", default=None)
@click.option('-c', '--confidence', help="Output confidence of mugration inference", default=False, is_flag=True)
@click.option('-m', '--missing_data', default="?", help="String indicating missing data")
@click.option('-o', '--out', help="Directory to write the output to")
@click.option('-v', '--verbose', default=1, type=int, help="Verbosity of output 0-6")
def mugration(name_column, states, attribute, tree, missing_data, rng_seed, frequencies, weights, power, verbose, confidence, out):
    outdir = get_outdir(out, '_wASR')
    
    tre = Phylo.read(tree, 'newick')
    taxon_in_tree = {n.name for n in tre.get_terminals()}
    
    # Parse states
    if not os.path.isfile(states):
        raise FileNotFoundError("States file does not exist")
    
    states_df = pd.read_csv(states, sep='\t' if states.lower().endswith('tsv') else ',',
                           skipinitialspace=True, na_filter=False)
    
    # Determine taxon name column
    if name_column:
        if name_column not in states_df.columns:
            raise ValueError(f"Column '{name_column}' not found in metadata file")
        taxon_name = name_column
    elif 'name' in states_df.columns:
        taxon_name = 'name'
    elif 'strain' in states_df.columns:
        taxon_name = 'strain'
    elif 'accession' in states_df.columns:
        taxon_name = 'accession'
    else:
        taxon_name = states_df.columns[0]
    
    print(f"Using column '{taxon_name}' as taxon name. This needs to match the taxa in the tree!")
    
    # Determine attribute column
    if attribute:
        if attribute not in states_df.columns:
            print(f"The specified attribute was not found in the metadata file {states}", file=sys.stderr)
            print(f"Available columns are: {', '.join(states_df.columns)}", file=sys.stderr)
            return 1
        attr = attribute
    else:
        attr = states_df.columns[1]
        print(f"Attribute for mugration inference was not specified. Using {attr}", file=sys.stderr)
    
    leaf_to_attr = {x[taxon_name]: str(x[attr]) if x[attr] else missing_data 
                    for _, x in states_df.iterrows() if x[taxon_name] in taxon_in_tree}
    
    weights_dict = parse_weights(leaf_to_attr, frequencies, weights, power)
    
    mug, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(
        tree, leaf_to_attr, missing_data=missing_data, verbose=verbose, 
        weights=weights_dict, rng_seed=rng_seed)
    
    unique_states = sorted(letter_to_state.values())
    
    # Output results
    print(f"\nCompleted mugration model inference of attribute '{attr}' for {tree}")
    
    basename = outdir
    gtr_name = basename + 'GTR.txt'
    with open(gtr_name, 'w', encoding='utf-8') as ofile:
        ofile.write('Character to attribute mapping:\n')
        for state in unique_states:
            ofile.write(f'  {reverse_alphabet[state]}: {state}\n')
        ofile.write(f'\n\n{str(mug.gtr)}\n')
    print(f"\nSaved inferred mugration model as: {gtr_name}")
    
    terminal_count = 0
    for n in mug.tree.find_clades():
        n.confidence = None
        if n.is_terminal() and len(n.name) > 40 and bioversion < "1.69":
            n.name = f"{n.name[:35]}_{terminal_count:03d}"
            terminal_count += 1
        n.comment = f'&{attr}="{letter_to_state[n.cseq[0]]}"'
    
    if confidence:
        conf_name = basename + 'confidence.csv'
        with open(conf_name, 'w', encoding='utf-8') as ofile:
            ofile.write('#name, ' + ', '.join(mug.gtr.alphabet) + '\n')
            for n in mug.tree.find_clades():
                ofile.write(n.name + ', ' + ', '.join([str(x) for x in n.marginal_profile[0]]) + '\n')
        print(f"Saved table with ancestral state confidences as: {conf_name}")
    
    outtree_name = basename + 'annotated_tree.nexus'
    Phylo.write(mug.tree, outtree_name, 'nexus')
    print(f"Saved annotated tree as: {outtree_name}")
    print("---Done!\n")
    
    return 0


if __name__ == '__main__':
    mugration()