# TransFact

TransFact is a Python script designed to estimate transmission patterns from pathogen transmission data. It uses a factorization approach to decompose transmission matrices into human mobility vectors and pathogen-specific source weights, optimized using a regularized objective function. The script processes transmission data from multiple files and outputs estimated human migration vectors and pathogen-specific features in CSV format.

## Prerequisites

To use TransFact, you need the following dependencies installed:

- Python 3.6+
- Required Python packages:
  - `numpy`
  - `pandas`
  - `scipy`
  - `click`

You can install the required packages using pip:

```bash
pip install numpy pandas scipy click
```

## Installation

1. Clone the repository:

```bash
git clone https://github.com/wASR_framework.git
cd wASR_framework/TransFact
```

2. Ensure all dependencies are installed as mentioned above.

3. Place the `TransFact.py` script in your working directory.

## Usage

The `TransFact.py` script processes one or more transmission files to estimate human mobility vectors and pathogen-specific source weights.

**Command-line options**:

- `-p, --prefix`: Prefix for output files (default: `transmission`).
- `-m, --minimum_involvement`: Minimum number of transmission events for a source or target to be included (default: 10).
- `transmission_files`: One or more input files containing transmission data.

**Example**:

```bash
python TransFact.py -p output -m 5 transmission1.txt transmission2.txt
```

**Output**:
- Human mobility vectors: `<prefix>.vectors.csv`
- Pathogen-specific source weights: `<prefix>.pathogens.csv`

## Input File Format

- **Transmission files**: Tab-separated text files where each line represents a transmission event with at least two columns:
  - First column: Source location (e.g., country or region).
  - Second column: Target location (e.g., country or region).
  - Lines with identical source and target locations are ignored.
  - Locations starting with "China" are normalized to "China".

**Example transmission file**:
```
USA	Canada
China	Japan
Canada	USA
```

## Workflow

1. **Parsing Transmission Data**:
   - Reads transmission files and filters sources and targets based on the `minimum_involvement` threshold.
   - Constructs transmission matrices for each input file, normalized by target sums.

2. **Matrix Factorization**:
   - Initializes human mobility vectors and pathogen-specific source weights.
   - Optimizes parameters using the L-BFGS-B method to minimize a loss function with regularization.
   - Applies a mask to exclude invalid pathogen-source combinations.

3. **Output**:
   - Saves human mobility vectors as a CSV file with targets as rows and sources as columns.
   - Saves pathogen-specific source weights as a CSV file with transmission files as rows and sources as columns.

## Notes

- The script uses a small epsilon value (`1e-10`) to avoid division by zero in normalization steps.
- The optimization process uses the `scipy.optimize.minimize` function with the L-BFGS-B method, configured for high precision (`ftol=1e-10`, `gtol=1e-8`).
- The `minimum_involvement` parameter ensures only frequently involved sources and targets are included in the analysis.
- The output CSV files are indexed by sorted target and source names for consistency.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue on GitHub for bug reports, feature requests, or improvements.
