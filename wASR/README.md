# wASR: Weighted Ancestral State Reconstruction

## Overview

The wASR package is designed for performing weighted ancestral state reconstruction (ASR) of discrete traits on phylogenetic trees. It provides a set of tools to analyze and infer the ancestral states of discrete characters, taking into account sampling biases and other factors.

## Installation

1. Clone the repository:
   ```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/wASR
   ```

2. Install Python dependencies:
   ```bash
   pip install ete3 click numpy pandas numba biopython
   ```


## Usage

### Command - Line Interface

The main entry point of the package is the `mugration` command in the `wASR.py` script, which can be used as follows:

Run the script:
```bash
python wASR.py -t <tree_file> -s <states_file> [options]
```

**Options**:
- `-t` / `--tree`: Required. Name of the file containing the tree in Newick, Nexus, or Phylip format. The branch length of the tree should be in units of the average number of nucleotide or protein substitutions per site. If no file is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed).
- `-s` / `--states`: Required. CSV or TSV file with discrete characters. The file should have a column for taxon names and a column for the discrete trait values.
- `-f` / `--frequencies`: CSV or TSV file with expected total counts of each state, such as the number of cases in each country.
- `-w` / `--weights`: CSV or TSV file with expected sampling rate of each state.
- `-p` / `--power`: Power adjustment of the weights. Default is 1.0.
- `-r` / `--rng_seed`: Random number generator seed for treetime.
- `-n` / `--name_column`: Label of the column to be used as the taxon name.
- `-a` / `--attribute`: Attribute to reconstruct. E.g., continent. Default is the second column of the states file.
- `-c` / `--confidence`: Output the confidence of mugration inference.
- `-m` / `--missing_data`: String indicating missing data. Default is ? or ''.
- `-o` / `--out`: Directory to write the output to.
- `-v` / `--verbose`: Verbosity of the output. Range is 0 - 6. Default is 1.


**How the command works**:
1. Tree and States Parsing:
  - The mugration function first reads the phylogenetic tree from the provided file using Bio.Phylo.read.
  - It then parses the states file, determining the taxon name column and the attribute column for reconstruction.

2. Weights Calculation:
  - The parse_weights function is called to calculate the weights for each taxon based on the provided frequencies or weights file, and the power adjustment.

3. Ancestral State Reconstruction:
  - The reconstruct_discrete_traits function is called to perform the actual ancestral state reconstruction. It sets up a General Time Reversible (GTR) model, initializes TreeAnc for the tree, and performs multiple iterations of inference and optimization.

4. Output Generation:
  - After the reconstruction, the inferred mugration model is saved to a GTR.txt file, including the character - to - attribute mapping.
  - If the --confidence option is provided, a confidence.csv file is generated with the ancestral state confidences for each node in the tree.
  - An annotated tree in Nexus format (annotated_tree.nexus) is saved, with the ancestral states annotated as node comments.


**Example**:
```bash
python wASR.py -t my_tree.newick -s my_states.csv -a continent -o output_dir
```

**Function - Level Usage**:
If you want to use the functions in your Python scripts, you can import them as follows:
```python
from wASR import reconstruct_discrete_traits, get_outdir, parse_weights

# Example usage of reconstruct_discrete_traits
import Bio.Phylo as Phylo
import pandas as pd
import numpy as np

tree = Phylo.read('my_tree.newick', 'newick')
states = pd.read_csv('my_states.csv')
leaf_to_attr = {x['name']: str(x['continent']) for xi, x in states.iterrows()}
weights = parse_weights(leaf_to_attr)

mug, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(tree, leaf_to_attr, weights=weights)
```


**Output**:
The mugration command generates the following output files:

1. GTR.txt: Contains the inferred mugration model and the mapping from characters to attributes.
2. confidence.csv: If the -c option is used, this file contains the ancestral state confidences for each node in the tree.
3. annotated_tree.nexus: A phylogenetic tree file with the ancestral states annotated as node comments.


**Configuration**:
The config.py file contains a set of global constants and parameters used throughout the project. These include:

- Numerical Precision: Constants like BIG_NUMBER, TINY_NUMBER, and SUPERTINY_NUMBER are used for numerical stability in calculations.
- Grid Sizes: Various grid sizes for numerical integration and distribution calculations, such as BRANCH_GRID_SIZE, NODE_GRID_SIZE, and N_INTEGRAL, are defined at different levels of precision (rough, normal, fine, ultra).
- Clock Tree Parameters: Parameters related to the clock tree, such as BRANCH_LEN_PENALTY, MAX_BRANCH_LENGTH, and coefficients for the autocorrelated molecular clock (MU_ALPHA and MU_BETA).

These parameters can be adjusted according to the specific requirements of your analysis.

## Contributing

Contributions are welcome! Please submit pull requests or open issues on GitHub for bug reports, feature requests, or improvements.
