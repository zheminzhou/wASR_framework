import os, sys
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio import __version__ as bioversion
from treeanc import TreeAnc
from gtr import GTR


def get_outdir(outdir, suffix='_wASR'):
    if outdir:
        if os.path.exists(outdir):
            if os.path.isdir(outdir):
                return outdir.rstrip('/') + '/'
            else:
                print("designated output location %s is not a directory"%outdir, file=sys.stderr)
        else:
            os.makedirs(outdir)
            return outdir.rstrip('/') + '/'

    from datetime import datetime
    outdir_stem = datetime.now().date().isoformat()
    outdir = outdir_stem + suffix.rstrip('/')+'/'
    count = 1
    while os.path.exists(outdir):
        outdir = outdir_stem + '-%04d'%count + suffix.rstrip('/')+'/'
        count += 1

    os.makedirs(outdir)
    return outdir

def parse_weights(leaf_to_attributes, frequencies=None, weights=None, power=0.5) :
    import collections
    if not weights :
        if frequencies :
            attribute_counts = collections.Counter(leaf_to_attributes.values())
            base_count = np.log(len(leaf_to_attributes.keys()))
            try :
                tmp_frequencies = pd.read_csv(frequencies, sep='\t' if frequencies[-3:]=='tsv' else ',', skipinitialspace=True)
                weight_dict = {}
                for state, freq in tmp_frequencies.values :
                    try :
                        freq = float(freq)
                        weight_dict[state] = freq/(base_count+attribute_counts.get(state, 0))
                    except :
                        pass
            except:
                raise ValueError(f"Loading of frequencies file '{frequencies}' failed!")
        else :
            weight_dict = {state:1. for state in leaf_to_attributes.values()}
    else :
        try:
            tmp_weights = pd.read_csv(weights, sep='\t' if weights[-3:]=='tsv' else ',', skipinitialspace=True)
            weight_dict = {}
            for state, weight in tmp_weights.values :
                try :
                    weight_dict[state] = float(weight)
                except :
                    pass
        except:
            raise ValueError("Loading of weights file '%s' failed!"%weights)
    weight_dict = { state:weight**power for state, weight in weight_dict.items() }
    mean_weight = np.mean([weight_dict[state] for state in leaf_to_attributes.values() if state in weight_dict])
    weights = {leaf:weight_dict.get(state, mean_weight)/mean_weight for leaf, state in leaf_to_attributes.items()}
    return weights


def reconstruct_discrete_traits(tree, traits, missing_data='?', pc=1.0, sampling_bias_correction=None,
                                weights=None, verbose=0, iterations=7, rng_seed=None):

    unique_states = sorted(set(traits.values()) - {missing_data})
    n_observed_states = len(unique_states)

    # make a map from states (excluding missing data) to characters in the alphabet
    # note that gap character '-' is chr(45) and will never be included here
    reverse_alphabet = {state:chr(65+i) for i, state in enumerate([state for state in unique_states if state!=missing_data])}
    alphabet = list(reverse_alphabet.values())
    # construct a look up from alphabet character to states
    letter_to_state = {v:k for k,v in reverse_alphabet.items()}


    # consistency checks
    if len(alphabet)<2:
        print("mugration: only one or zero states found -- this doesn't make any sense", file=sys.stderr)
        return None, None, None

    n_states = len(alphabet)
    missing_char = chr(65+n_states)
    reverse_alphabet[missing_data]=missing_char
    letter_to_state[missing_char]=missing_data

    ###########################################################################
    ### construct gtr model
    ###########################################################################

    # set up dummy matrix
    W = np.ones((n_states,n_states), dtype=float)

    mugration_GTR = GTR.custom(pi = None, W=W, alphabet = np.array(alphabet))
    mugration_GTR.profile_map[missing_char] = np.ones(n_states)
    mugration_GTR.ambiguous=missing_char

    ###########################################################################
    ### set up treeanc
    ###########################################################################
    treeanc = TreeAnc(tree, gtr=mugration_GTR, verbose=verbose, ref='A',
                      convert_upper=False, one_mutation=0.0001, rng_seed=rng_seed)
    treeanc.use_mutation_length = False
    pseudo_seqs = {n.name: {0:reverse_alphabet[traits[n.name]] if n.name in traits else missing_char}
                   for n in treeanc.tree.get_terminals()}

    # construct the vector with weights to be used as equilibrium frequency

    valid_seq = np.array([s[0]!=missing_char for s in pseudo_seqs.values()])
    print("Assigned discrete traits to %d out of %d taxa.\n"%(np.sum(valid_seq),len(valid_seq)))
    treeanc.aln = pseudo_seqs
    # treeanc.weights = weights
    for node in treeanc.tree.get_terminals() :
        node.mask = weights[node.name] if weights else 1.
    # for node in treeanc.tree.get_nonterminals() :
    #     node.mask = 1.

    try:
        ndiff = treeanc.infer_ancestral_sequences(method='ml', infer_gtr=True,
            store_compressed=False, pc=pc, marginal=True, normalized_rate=False,
            fixed_pi=None, reconstruct_tip_states=True)
        treeanc.optimize_gtr_rate()
    except :
        raise "\nAncestral reconstruction failed, please see above for error messages and/or rerun with --verbose 4\n"

    for i in range(iterations):
        treeanc.infer_gtr(marginal=True, normalized_rate=False, pc=pc, fixed_pi=None)
        treeanc.optimize_gtr_rate()

    if sampling_bias_correction:
        treeanc.gtr.mu *= sampling_bias_correction

    ndiff = treeanc.infer_ancestral_sequences(infer_gtr=False, store_compressed=False,
                                 marginal=True, normalized_rate=False,
                                 reconstruct_tip_states=True)

    return treeanc, letter_to_state, reverse_alphabet


import click

@click.command()
@click.option('-t', '--tree', required = True, type=str, help="Name of file containing the tree in newick, nexus, or phylip format, the branch length of the tree should be in units of average number of nucleotide or protein substitutions per site. If no file is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed).")
@click.option('-s', '--states', required = True, type=str, help="csv or tsv file with discrete characters. #name,continent,country taxon1,oceania,micronesia ...")
@click.option('-f', '--frequencies', help="csv or tsv file with expected total counts of each state, such as number of cases in each country. E.g.: #country,weight\nmicronesia,0.1 ...")
@click.option('-w', '--weights', help="csv or tsv file with expected sampling rate of each state. E.g.: #country,weight\nmicronesia,0.1 ...")
@click.option('-p', '--power', help="power adjustment of the weights default:1.0", default=1.0)
@click.option('-r', '--rng_seed', type=int, help="random number generator seed for treetime")
@click.option('-n', '--name_column', help="label of the column to be used as taxon name")
@click.option('-a', '--attribute', help="attribute to reconstruct. e.g. continent", default=None)
@click.option('-c', '--confidence', help="output confidence of mugration inference", default=False, is_flag=True)
@click.option('-m', '--missing_data', default="?", help="string indicating missing data. default:? or ''")
@click.option('-o', '--out', help="directory to write the output to")
@click.option('-v', '--verbose', default=1, type=int, help="verbosity of output 0-6")
def mugration(name_column, states, attribute, tree, missing_data, rng_seed, frequencies, weights, power, verbose, confidence, out):
    outdir = get_outdir(out, '_wASR')

    tre = Phylo.read(tree, 'newick')
    taxon_in_tree = set([n.name for n in tre.get_terminals()])

    ###########################################################################
    ### Parse states
    ###########################################################################
    if os.path.isfile(states):
        states = pd.read_csv(states, sep='\t' if states.lower().endswith('tsv') else ',',
                             skipinitialspace=True, na_filter=False)
    else:
        raise ("file with states does not exist")

    if name_column:
        if name_column in states.columns:
            taxon_name = name_column
        else:
            raise ("Error: specified column '%s' for taxon name not found in meta data file with columns: "%name_column + " ".join(states.columns))

    elif 'name' in states.columns: taxon_name = 'name'
    elif 'strain' in states.columns: taxon_name = 'strain'
    elif 'accession' in states.columns: taxon_name = 'accession'
    else:
        taxon_name = states.columns[0]
    print("Using column '%s' as taxon name. This needs to match the taxa in the tree!"%taxon_name)

    if attribute:
        if attribute in states.columns:
            attr = attribute
        else:
            print("The specified attribute was not found in the metadata file "+states, file=sys.stderr)
            print("Available columns are: "+", ".join(states.columns), file=sys.stderr)
            return 1
    else:
        attr = states.columns[1]
        print("Attribute for mugration inference was not specified. Using " + attr, file=sys.stderr)

    leaf_to_attr = {x[taxon_name]:str(x[attr]) if x[attr] else missing_data for xi, x in states.iterrows() if x[taxon_name] in taxon_in_tree}

    weights = parse_weights(leaf_to_attr, frequencies, weights, power)

    mug, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(tree, leaf_to_attr,
                missing_data=missing_data, verbose=verbose, weights=weights, rng_seed=rng_seed)

    unique_states = sorted(letter_to_state.values())
    ###########################################################################
    ### output
    ###########################################################################
    print("\nCompleted mugration model inference of attribute '%s' for"% attr, tree)

    basename = outdir
    gtr_name = basename + 'GTR.txt'
    with open(gtr_name, 'w', encoding='utf-8') as ofile:
        ofile.write('Character to attribute mapping:\n')
        for state in unique_states:
            ofile.write('  %s: %s\n'%(reverse_alphabet[state], state))
        ofile.write('\n\n'+str(mug.gtr)+'\n')
        print("\nSaved inferred mugration model as:", gtr_name)

    terminal_count = 0
    for n in mug.tree.find_clades():
        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        n.comment= '&%s="'%attr + letter_to_state[n.cseq[0]] +'"'

    if confidence:
        conf_name = basename+'confidence.csv'
        with open(conf_name, 'w', encoding='utf-8') as ofile:
            ofile.write('#name, '+', '.join(mug.gtr.alphabet)+'\n')
            for n in mug.tree.find_clades():
                ofile.write(n.name + ', '+', '.join([str(x) for x in n.marginal_profile[0]])+'\n')
        print("Saved table with ancestral state confidences as:", conf_name)

    # write tree to file
    outtree_name = basename+'annotated_tree.nexus'
    Phylo.write(mug.tree, outtree_name, 'nexus')
    print("Saved annotated tree as:", outtree_name)
    print("---Done!\n")

    return 0



if __name__ == '__main__' :
    mugration()
