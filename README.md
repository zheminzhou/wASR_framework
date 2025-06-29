# wASR Framework: A Comprehensive Toolkit for Genomic and Phylogenetic Analysis

## Overview
The wASR framework is a collection of Python and R scripts designed to perform a wide range of genomic and phylogenetic analyses. It includes tools for ancestral state reconstruction, trait simulation, mutation analysis, and more. Each tool in the framework is modular and can be used independently or in combination with other tools to build complex analysis pipelines.

## Tools in the Framework

### 1. wASR: Weighted Ancestral State Reconstruction
- **Purpose**: Reconstruct ancestral states with weighted probabilities, considering the sampling rate and frequency of each state.
- **Installation**:
```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/wASR
pip install ete3 click numpy pandas numba biopython
```

### 2. ASR_simulations: Trait Simulation and Ancestral State Reconstruction Tools
- **Purpose**: Simulate traits on phylogenetic trees and perform ancestral state reconstruction using various methods.
- **Installation**:
```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/ASR_simulations
pip install ete3 click numpy pandas numba
```
```R
install.packages("ape")
```
Follow the instructions for [TreeTime](https://github.com/neherlab/treetime#installation) and download and set up [wASR](https://github.com/EPHI-bioinformatics/wASR) as per its documentation.


### 3. DHMM: Divergent Hidden Markov Model
- **Purpose**: Identify diversified regions in genomic sequences using a Hidden Markov Model (HMM) approach.
- **Installation**:
```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/DHMM
pip install numpy pandas numba
```
- **Citation**: Z. Zhou, A. McCann, F. Weill, C. Blin, S. Nair, J. Wain, G. Dougan, & M. Achtman, *Transient Darwinian selection in Salmonella enterica serovar Paratyphi A during 450 years of global spread of enteric fever*, Proc. Natl. Acad. Sci. U.S.A. 111 (33) 12199-12204, [https://doi.org/10.1073/pnas.1411012111](https://doi.org/10.1073/pnas.1411012111) (2014).

### 4. TransFact
- **Purpose**: Analyze transmission factors in genomic data.
- **Installation**:
```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/TransFact
pip install numpy pandas scipy click
```

### 5. NNPred
- **Purpose**: Nearest Neighbor-based trait source prediction tool.
- **Installation**:
```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/NNPred
```
Ensure all dependencies are installed as mentioned above. Place the `NNPred_train.py` and `NNPred_query.py` scripts in your working directory.


### 6. phylogeny_transmission: Phylogenetic Trait Subsampling and Transmission Analysis
- **Purpose**: Perform phylogenetic trait subsampling and transmission analysis.
- **Installation**:
```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/phylogeny_transmission
pip install ete3 click numpy pandas numba
```
Ensure the `ete3_extensions.py` module (from the previous package) is available in the same directory or Python path, as it is required for Nexus file handling.

### 7. mutation_analysis
- **Purpose**: Analyze mutation hotspots in genomic data, associate mutations with genes and geographic traits, perform enrichment analyses, and visualize the results.
- **Dependencies**:
```bash
pip install click pandas numpy scipy ete3 pycirclize matplotlib keggtools statsmodels
```
- **Pipeline Workflow**:
    1. Annotate Diversified Regions with Genes: Run `DivHMM_region_annotation.py` to associate diversified regions with genes and mutation counts.
    2. Add KO Annotations: Run `DivHMM_to_genes.py` to extend annotations with KO identifiers and categorize mutations by region type.
    3. Perform KO Enrichment Analysis: Run `DivHMM_to_genes_KO_enrichment.py` to identify enriched KO categories in diversified and recombining regions.
    4. Associate Mutations with Traits: Run `DivHMM_to_trait.py` to associate mutations with geographic or other traits using phylogenetic data.
    5. Map Genes to Enriched Geographic Regions: Run `DivHMM_to_trait_gene_enrichment.py` to link genes and KO terms to geographically enriched regions.
    6. Analyze Temporal Mutation Rates: Run `DivHMM_to_timed.py` to analyze mutation rates over time (optional, if temporal analysis is needed).
    7. Visualize Results: Run `DivHMM_plot.py` to generate a circular genome plot visualizing genes, mutation densities, and geographic associations.
    8. Compute Geographic Region Similarities: Run `DivHMM_to_trait_region_network.py` to calculate similarities between geographic regions based on shared mutation hotspots.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue on GitHub for bug reports, feature requests, or improvements.

## License
This project is licensed under the GPL-3.0 License. See the respective `LICENSE` files in each subdirectory for details.

## Contact
For questions or support regarding the entire framework or specific tools, please open an issue on the GitHub repository or contact the maintainers.
