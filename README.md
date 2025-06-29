# wASR Framework: A Modular Toolkit for Genomic and Phylogenetic Analysis

## Overview
The wASR Framework is a modular collection of Python and R tools designed for genomic and phylogenetic analyses, including ancestral state reconstruction, trait simulation, mutation analysis, and transmission factor analysis. Each tool is independent, allowing flexible use in standalone or integrated workflows.

## Installation
Clone the repository to get started:
```bash
git clone https://github.com/zheminzhou/wASR_framework.git
cd wASR_framework
```

Install core dependencies for Python-based tools:
```bash
pip install ete3 click numpy pandas numba biopython scipy matplotlib keggtools statsmodels
```

For R-based tools, install the required package:
```R
install.packages("ape")
```

Some tools may require additional setup (e.g., [TreeTime](https://github.com/neherlab/treetime#installation)). Refer to specific tool sections for details.

## Tools

### 1. wASR: Weighted Ancestral State Reconstruction
- **Purpose**: Reconstructs ancestral states with weighted probabilities based on sampling rates and state frequencies.
- **Setup**:
  ```bash
  cd wASR
  ```
  Dependencies are covered by the core installation.

### 2. ASR_simulations: Trait Simulation and Reconstruction
- **Purpose**: Simulates traits on phylogenetic trees and performs ancestral state reconstruction.
- **Setup**:
  ```bash
  cd ASR_simulations
  ```
  Follow [TreeTime](https://github.com/neherlab/treetime#installation) and [wASR](https://github.com/EPHI-bioinformatics/wASR) documentation for additional setup.

### 3. DHMM: Divergent Hidden Markov Model
- **Purpose**: Identifies diversified regions in genomic sequences using a Hidden Markov Model.
- **Setup**:
  ```bash
  cd DHMM
  ```
- **Citation**: Z. Zhou et al., *Transient Darwinian selection in Salmonella enterica serovar Paratyphi A during 450 years of global spread of enteric fever*, Proc. Natl. Acad. Sci. U.S.A. 111(33), 12199-12204, 2014. [DOI:10.1073/pnas.1411012111](https://doi.org/10.1073/pnas.1411012111).

### 4. TransFact: Transmission Factor Analysis
- **Purpose**: Analyzes transmission factors in genomic data.
- **Setup**:
  ```bash
  cd TransFact
  ```

### 5. NNPred: Nearest Neighbor Trait Prediction
- **Purpose**: Predicts trait sources using a nearest neighbor approach.
- **Setup**:
  ```bash
  cd NNPred
  ```
  Ensure `NNPred_train.py` and `NNPred_query.py` are in your working directory.

### 6. phylogeny_transmission: Phylogenetic Trait Subsampling
- **Purpose**: Conducts phylogenetic trait subsampling and transmission analysis.
- **Setup**:
  ```bash
  cd phylogeny_transmission
  ```
  Ensure `ete3_extensions.py` is in the Python path for Nexus file handling.

### 7. mutation_analysis: Mutation Hotspot Analysis
- **Purpose**: Analyzes mutation hotspots, associates mutations with genes and traits, and visualizes results.
- **Workflow**:
  1. **Annotate Diversified Regions**: Run `DivHMM_region_annotation.py` to link regions with genes and mutation counts.
  2. **Add KO Annotations**: Run `DivHMM_to_genes.py` to include KO identifiers and categorize mutations.
  3. **KO Enrichment Analysis**: Run `DivHMM_to_genes_KO_enrichment.py` to identify enriched KO categories.
  4. **Trait Association**: Run `DivHMM_to_trait.py` to associate mutations with geographic or other traits.
  5. **Geographic Enrichment**: Run `DivHMM_to_trait_gene_enrichment.py` to link genes/KO terms to enriched regions.
  6. **Temporal Analysis** (Optional): Run `DivHMM_to_timed.py` for mutation rate analysis over time.
  7. **Region Similarity**: Run `DivHMM_to_trait_region_network.py` to compute similarities between geographic regions.

## Contributing
We welcome contributions! Submit pull requests or open issues on [GitHub](https://github.com/zheminzhou/wASR_framework) for bugs, features, or improvements.

## License
Licensed under the GPL-3.0 License. See `LICENSE` files in each subdirectory for details.

## Contact
For support or inquiries, open an issue on the [GitHub repository](https://github.com/zheminzhou/wASR_framework) or contact the maintainers.
