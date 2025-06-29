# Phylogenetic Trait Subsampling and Transmission Analysis

This package provides tools for subsampling phylogenetic trees based on trait distributions, summarizing ancestral state reconstructions across subsampled datasets, and identifying potential transmission events in phylogenetic trees. It is designed for researchers working with phylogenetic data to analyze trait distributions and infer transmission patterns.

## Overview

The package includes three main scripts:
- **`EPHI_trait_subsample.py`**: Subsamples a phylogenetic tree based on trait data, ensuring balanced or weighted representation of traits.
- **`EPHI_trait_subsample_summary.py`**: Aggregates ancestral state reconstruction results across multiple subsampled trees to determine consensus states and their probabilities.
- **`tree2transmission.py`**: Identifies potential transmission events in a tree by comparing ancestral states between nodes, with configurable criteria for state confidence and node inclusion.
- **`ete3_extensions.py`**: A utility module (assumed to be included from a previous context) providing functions for reading, writing, and manipulating phylogenetic trees in Nexus format.

These tools leverage Python and the `ete3` library for tree manipulation, with support for handling large datasets and complex trait distributions.

## Prerequisites

To use this package, you need the following software and libraries installed:

- **Python 3.x** with:
  - `ete3`: For phylogenetic tree manipulation.
  - `click`: For command-line interface creation.
  - `numpy`: For numerical operations.
  - `pandas`: For data handling.
  - `numba`: For optimized numerical computations (used in `ete3_extensions.py`).
- **Conda** (recommended) for managing Python environments.
- A Unix-like environment (Linux or macOS) is recommended, as the scripts use standard file operations.

## Installation

1. Clone the repository:
   ```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/phylogeny_transmission
   ```

2. Install Python dependencies:
   ```bash
   pip install ete3 click numpy pandas numba
   ```

3. (Optional) Set up a Conda environment:
   ```bash
   conda create -n phylo_analysis python=3.8
   conda activate phylo_analysis
   conda install -c conda-forge ete3 click numpy pandas numba
   ```

4. Ensure the `ete3_extensions.py` module (from the previous package) is available in the same directory or Python path, as it is required for Nexus file handling.

## Usage

### 1. Subsample Phylogenetic Trees by Traits

The `EPHI_trait_subsample.py` script subsamples a phylogenetic tree based on trait data, either by limiting the number of samples per trait or using weights to determine sampling probabilities.

Run the script:
```bash
python EPHI_trait_subsample.py -d <output_directory> -n <nexus_file> -t <trait_file> [-w <weight>]
```

**Options**:
- `-d` / `--outdir`: Output directory for results.
- `-n` / `--nexus`: Input Nexus file containing the phylogenetic tree.
- `-t` / `--trait`: Tab-separated file with trait data (columns: node ID, trait).
- `-w` / `--weight`: Either an integer (maximum number of samples per trait, default: 10) or a tab-separated file with weights (columns: trait, weight).

**Example**:
```bash
python EPHI_trait_subsample.py -d output -n tree.nex -t traits.txt -w 5
```

**Output**:
- A directory `<output_directory>` containing:
  - `dating.out.nex`: The original tree in Nexus format.
  - `trait.txt`: A tab-separated file with subsampled node IDs and their traits.

### 2. Summarize Ancestral State Reconstructions

The `EPHI_trait_subsample_summary.py` script aggregates ancestral state reconstruction results from multiple subsampled trees, computing consensus states and their probabilities for each node.

Run the script:
```bash
python EPHI_trait_subsample_summary.py -o <output_nexus> -l <list_file>
```

**Options**:
- `-o` / `--out_nexus`: Output Nexus file with summarized annotations.
- `-l` / `--list`: Text file listing paths to subsampled Nexus files.

**Example**:
```bash
python EPHI_trait_subsample_summary.py -o summary.nex -l subsamples.txt
```

**Input**:
- `subsamples.txt`: A text file with one Nexus file path per line, each containing a tree with `traits.prop` and `traits.list` annotations.

**Output**:
- A Nexus file (`summary.nex`) with a tree where each node has annotations:
  - `state`: The most likely trait.
  - `state.prop`: The normalized probability of the most likely trait.
  - `traits.list`: List of all traits observed for the node.
  - `traits.prop`: Corresponding probabilities for each trait.

### 3. Identify Transmission Events

The `tree2transmission.py` script analyzes a tree to identify potential transmission events by comparing ancestral states between parent and child nodes, filtering based on state confidence.

Run the script:
```bash
python tree2transmission.py [-p <proportion>] [-a] <nexus_file>
```

**Options**:
- `-p` / `--proportion`: Minimum probability threshold for a node's state to be considered (default: 0.5).
- `-a` / `--all_nodes`: If set, include all nodes in the output, not just those involved in transmissions.
- `<nexus_file>`: Input Nexus file with a tree containing `state`, `state.prop`, and `date` annotations.

**Example**:
```bash
python tree2transmission.py -p 0.7 -a summary.nex
```

**Output**:
- Tab-separated output to stdout with columns:
  - Parent state
  - Child state
  - Parent state probability
  - Child state probability
  - Parent date
  - Child date
  - Number of descendants (child node)
  - Parent node name
  - Child node name

### 4. Utility Functions

The `ete3_extensions.py` module (assumed from the previous package) provides essential functions:
- `read_nexus`: Reads trees and annotations from Nexus files.
- `write_nexus`: Writes trees with annotations to Nexus files.
- `prune`: Prunes trees to retain specific tips.
- `pairwise_distance`: Computes pairwise distances between leaves.

Ensure this module is available in your working directory or Python path.

## Directory Structure

After running the scripts, the directory structure may look like this:

```
phylogenetic-trait-analysis/
├── output/
│   ├── dating.out.nex
│   ├── trait.txt
│   └── ...
├── summary.nex
├── subsamples.txt
├── EPHI_trait_subsample.py
├── EPHI_trait_subsample_summary.py
├── tree2transmission.py
├── ete3_extensions.py
└── README.md
```

## Notes

- The scripts assume that input Nexus files contain trees with appropriate annotations (e.g., `traits.prop`, `traits.list`, `state`, `state.prop`, `date`).
- The `EPHI_trait_subsample.py` script supports two modes for subsampling: fixed sample size per trait or weighted sampling based on a provided weight file.
- The `tree2transmission.py` script requires `date` annotations, which may need to be added by other tools (e.g., TreeTime for dating).
- Temporary directories created by `EPHI_trait_subsample.py` are not automatically deleted. Use `rm -rf <outdir>` to clean up if needed.
- Ensure sufficient disk space for large trees or multiple subsamples.

## Contributing

Contributions are welcome! Please submit pull requests or open issues on GitHub for bug reports, feature requests, or improvements.
