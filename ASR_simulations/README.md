# Trait Simulation and Ancestral State Reconstruction Tools

This package provides a suite of tools for simulating trait transitions on phylogenetic trees and evaluating ancestral state reconstruction (ASR) methods. It generates simulated datasets, subsamples trees, and compares the performance of various ASR tools, including TreeTime, wASR, and APE (with and without weighting).

## Overview

The package includes scripts to:
- Generate random phylogenetic trees and simulate discrete trait transitions with varying frequencies using R (`trait_sim.R`).
- Subsample trees to create datasets with specific trait frequency ratios (`subsample_sim.py`).
- Compare the performance of ASR tools (TreeTime, wASR, APE, and weighted APE) on simulated data (`compare_tools.py`).
- Evaluate the impact of incomplete or inaccurate weights in weighted ASR (`wASR_incomplete_weights.py`).
- Provide utility functions for tree manipulation and nexus file handling (`ete3_extensions.py`).
- Implement a weighted ancestral character estimation method (`wACE.R`).
- Summarize and analyze results across different simulation scenarios (`compare_tools.summary.py`).

The tools are designed to work with Python, R, and external dependencies like TreeTime and wASR, and they produce output in formats suitable for further analysis (e.g., Nexus, CSV).

## Prerequisites

To use this package, you need the following software and libraries installed:

- **Python 3.x** with:
  - `ete3`: For phylogenetic tree manipulation.
  - `click`: For command-line interface creation.
  - `numpy`: For numerical operations.
  - `pandas`: For data handling.
  - `numba`: For optimized numerical computations.
- **R** with:
  - `ape`: For phylogenetic analysis and ancestral state reconstruction.
- **External Tools**:
  - [TreeTime](https://github.com/neherlab/treetime): For ancestral state reconstruction.
  - [wASR](https://github.com/EPHI-bioinformatics/wASR): Weighted ancestral state reconstruction tool.
- **Conda** (recommended) for managing Python and R environments.
- A Unix-like environment (Linux or macOS) is assumed for running scripts, as some paths are hardcoded (e.g., `/home/zhemin/miniconda3/envs/p312/bin/treetime`).

## Installation

1. Clone the repository:
   ```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/ASR_simulations
   ```

2. Install Python dependencies:
   ```bash
   pip install ete3 click numpy pandas numba
   ```

3. Install R dependencies:
   ```R
   install.packages("ape")
   ```

4. Install TreeTime and wASR:
   - Follow the instructions for [TreeTime](https://github.com/neherlab/treetime#installation).
   - wASR is distributed as part of this package.

5. (Optional) Set up a Conda environment for Python and R:
   ```bash
   conda create -n trait_sim python=3.8 r-base
   conda activate trait_sim
   conda install -c conda-forge ete3 click numpy pandas numba r-ape
   ```

6. Update hardcoded paths in `compare_tools.py`, `wASR_incomplete_weights.py`, and `compare_tools.summary.py` to match your local installations of TreeTime and wASR.

## Usage

### 1. Simulate Phylogenetic Trees and Traits

The `trait_sim.R` script generates random phylogenetic trees with simulated discrete traits under different frequency scenarios (e.g., 1:1, 3:1, 10:1, 30:1, 100:1). Each tree has 10,000 tips, and multiple simulations are run per tree.

Run the script:
```bash
Rscript trait_sim.R
```

**Output**:
- Nexus files in the `phylo_simulations0` directory, named as `tree_<tree_id>_<scenario>_sim_<sim_id>.nex`.
- Each file contains a tree and trait assignments for tips and internal nodes.

### 2. Subsample Trees

The `subsample_sim.py` script subsamples the simulated trees to create datasets with specific trait frequency ratios (e.g., 1:1, 1:3, 1:10, 1:30) and a maximum of 200 tips per state.

Run the script:
```bash
python subsample_sim.py <nexus_file>
```

**Example**:
```bash
python subsample_sim.py phylo_simulations0/tree_1_x1_1_sim_1.nex
```

**Output**:
- Subsampled Nexus files named `<prefix>.sub_<p0>_<p1>.ite_<ite>.nex`, where `<prefix>` is the input file prefix, `<p0>_<p1>` are the frequency ratios, and `<ite>` is the iteration number.

### 3. Compare ASR Tools

The `compare_tools.py` script evaluates the performance of TreeTime, wASR, APE (SYM model), and weighted APE on subsampled trees. It compares reconstructed ancestral states against ground truth and reports accuracy metrics.

Run the script:
```bash
python compare_tools.py -n <nexus_file> -t <tasks>
```

**Options**:
- `-n`: Path to the subsampled Nexus file.
- `-t`: Tasks to run (e.g., `all`, `treetime`, `wASR`, `ape`, `wAPE`).

**Example**:
```bash
python compare_tools.py -n phylo_simulations0/tree_1_x1_1_sim_1.sub_1_1.ite_0.nex -t all
```

**Output**:
- A directory named after the input Nexus file prefix containing:
  - `tree.nwk`: The subsampled tree in Newick format.
  - `states.csv`: Trait assignments for tips.
  - `states.weight`: Weights for wASR and weighted APE (based on frequency ratios).
  - Output files from each tool (e.g., `treetime.out`, `wASR.out`, `SYM_result.csv`).
- Printed results showing accuracy metrics for each tool.

### 4. Evaluate Incomplete Weights

The `wASR_incomplete_weights.py` script tests the robustness of wASR to inaccurate weights by varying the weight of one state across a range of error values.

Run the script:
```bash
python wASR_incomplete_weights.py -n <nexus_file>
```

**Example**:
```bash
python wASR_incomplete_weights.py -n phylo_simulations0/tree_1_x1_1_sim_1.sub_1_1.ite_0.nex
```

**Output**:
- A directory named after the input Nexus file prefix containing:
  - `tree.nwk`, `states.csv`, and `states.weight` (as above).
  - wASR output files for each weight variation (`iWeight_<err>.out`).
- Printed results showing accuracy metrics for each weight variation.

### 5. Summarize Results

The `compare_tools.summary.py` script aggregates results from `compare_tools.py` across multiple simulations and reports mean and standard deviation of accuracy for each tool and scenario.

Run the script:
```bash
python compare_tools.summary.py -r <result_file>
```

**Example**:
```bash
python compare_tools.summary.py -r results.txt
```

**Output**:
- Tab-separated output with columns: scenario (mutation rate, ancestor frequency, subsampling frequency), tool, mean accuracy, and standard deviation.

### 6. Utility Functions

The `ete3_extensions.py` module provides functions for:
- Computing pairwise distances between leaves (`pairwise_distance`).
- Pruning trees to retain specific tips (`prune`).
- Reading and writing Nexus files with annotations (`read_nexus`, `write_nexus`).
- Iterating over trees in various formats (`iter_trees`, `read_trees`).

### 7. Weighted Ancestral Character Estimation

The `wACE.R` script implements a weighted version of the `ace` function from the `ape` package, allowing weights to be applied to tip states during ancestral state reconstruction.

This is used internally by `compare_tools.py` for the `wAPE` task.

## Directory Structure

After running the scripts, the directory structure will look like this:

```
trait-simulation-asr/
├── phylo_simulations0/
│   ├── tree_1_x1_1_sim_1.nex
│   ├── tree_1_x1_1_sim_1.sub_1_1.ite_0.nex
│   ├── tree_1_x1_1_sim_1.sub_1_1.ite_0/
│   │   ├── tree.nwk
│   │   ├── states.csv
│   │   ├── states.weight
│   │   ├── treetime.out/
│   │   ├── wASR.out/
│   │   ├── SYM_result.csv
│   │   ├── ape.script
│   │   ├── ape_weighted.script
│   │   └── ...
├── trait_sim.R
├── subsample_sim.py
├── compare_tools.py
├── compare_tools.summary.py
├── wASR_incomplete_weights.py
├── ete3_extensions.py
├── wACE.R
└── README.md
```

## Notes

- The scripts assume specific file paths for TreeTime and wASR. Update these in the Python scripts if your installation paths differ.
- The `trait_sim.R` script generates a large number of files (100 trees × 5 scenarios × 10 simulations = 5000 Nexus files). Ensure sufficient disk space.
- The `compare_tools.py` and `wASR_incomplete_weights.py` scripts create temporary directories for each input file, which are not automatically deleted (though `compare_tools.py` has a commented-out `shutil.rmtree` call for cleanup).
- The weighted APE implementation (`wACE.R`) assumes a symmetric (SYM) model for discrete traits. The commented-out ARD model code can be enabled if needed.
- Results are reported as counts of correct/incorrect predictions with high/low confidence, based on a probability threshold of 0.7.

## Contributing

Contributions are welcome! Please submit pull requests or open issues on GitHub for bug reports, feature requests, or improvements.
