# int-iScore

**Int**egrated **iScore** - A comprehensive toolkit for protein complex interface analysis from AI predictions and molecular dynamics.

## Overview

int-iScore provides a unified framework for evaluating protein-protein interaction quality from various AI prediction outputs (AlphaFold2, AlphaFold3) and molecular dynamics trajectories. It calculates multiple interface quality metrics including DockQ, SOAP, DOPE, CPscore, buried surface area, clashes, and Ramachandran analysis.

## Features

- **Multi-mode Analysis**: Support for AlphaFold3 with/without reference and MD trajectory snapshots
- **RMSD-based Clustering**: Automatic selection of representative structures from large datasets
- **MPI Parallelization**: Efficient processing of large model collections
- **Comprehensive Metrics**: DockQ, SOAP, DOPE, CPscore, pLDDT, BSA, clashes, Ramachandran
- **CIF/PDB Support**: Handles both AlphaFold3 CIF outputs and traditional PDB files

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/lluoto/int-iScore.git
cd int-iScore

# Create conda environment (Python 3.10+ required for freesasa)
conda env create -f env.yml
conda activate int_iscore

# Install the package
pip install -e .

# Extract required SVM models
tar -jxvf cp_svm.tar.bz2
tar -jxvf frustratometer.tar.bz2

# Shape Complementarity is calculated using a standalone Python implementation
# and does not require any additional external dependencies
```

### Dependencies

**Core Requirements:**
- Python 3.10+ (required for freesasa)
- numpy >= 1.20.0
- scipy >= 1.7.0
- biopython >= 1.79
- freesasa >= 2.1.0
- tqdm >= 4.62.0
- matplotlib >= 3.5.0
- scikit-learn >= 1.0.0

**Optional:**
- mpi4py >= 3.1.0 (for parallel processing)
- modeller >= 10.0 (for SOAP/DOPE scoring)

**External Software:**
- Qhull (for Voronoi tessellation in intercaat)
- DockQ (included in package)
- Modeller license (for SOAP scoring)

## Quick Start

### Mode 1: AlphaFold3 with Reference (DockQ calculation)

```bash
# Binary interface (A interacts with B)
int-iscore af3-ref \
    -i /path/to/af3/output \
    -c A B \
    -n my_analysis \
    -o results/

# Multi-chain interface (A interacts with B and C) - using first chain as query
int-iscore af3-ref \
    -i /path/to/af3/output \
    -c A B C \
    -n my_analysis \
    -o results/

# Complex interface (A interacts with B, C, D, E) - using first chain as query
int-iscore af3-ref \
    -i /path/to/af3/output \
    -c A B C D E \
    -n my_analysis \
    -o results/

# Multiple query chains (A and B both interact with C, D, E) - use ':' separator
int-iscore af3-ref \
    -i /path/to/af3/output \
    -c A B : C D E \
    -n my_analysis \
    -o results/
```

### Mode 2: Molecular Dynamics Snapshots

```bash
# Binary interface
int-iscore md \
    -i /path/to/trajectories \
    -c A B \
    --cluster-method spectral \
    --n-clusters 20 \
    -n md_analysis

# Multi-chain interface (A interacts with B and C)
int-iscore md \
    -i /path/to/trajectories \
    -c A B C \
    --cluster-method spectral

# Multiple query chains (A,B both interact with C,D,E)
int-iscore md \
    -i /path/to/trajectories \
    -c A B : C D E \
    --cluster-method spectral
```

### Mode 3: AlphaFold3 without Reference

```bash
# Binary interface
int-iscore af3-nonref \
    -i /path/to/af3/output \
    -c A B \
    --format cif \
    -n my_analysis

# Multi-chain interface (A interacts with B, C, D)
int-iscore af3-nonref \
    -i /path/to/af3/output \
    -c A B C D \
    --format cif

# Multiple query chains (A,B both interact with C,D,E)
int-iscore af3-nonref \
    -i /path/to/af3/output \
    -c A B : C D E \
    --format cif
```

## Project Structure

```
int-iScore/
├── src/int_iscore/
│   ├── __init__.py
│   ├── cli.py                 # Command-line interface
│   ├── core/
│   │   ├── __init__.py
│   │   ├── metrics.py         # Core metric functions
│   │   └── parser.py         # Argument parsing
│   ├── modes/
│   │   ├── __init__.py
│   │   ├── af3_reference.py   # AF3 with reference mode
│   │   ├── af3_nonref.py     # AF3 without reference mode
│   │   └── md_snapshots.py   # MD trajectory mode
│   └── utils/
│       ├── __init__.py
│       ├── sc_ec.py           # SC calculation interface
│       └── sc_calculator.py  # Standalone SC implementation
├── intercaat_functions.py     # Intercaat contact detection
├── intercaat_config.ini       # Intercaat configuration
├── svmlight/                  # CPscore SVM models
├── reference/                 # Reference structures
└── README.md
```

## Interface Metrics

| Feature | Range | Description |
|---------|-------|-------------|
| DockQ | [0, 1] | Global docking quality vs reference |
| SOAP | [0, 1] | Statistical potential score |
| DOPE | [0, 1] | Discrete optimized protein energy |
| CPscore | [0, 1] | Contact preference score |
| BSA | [0, ∞] | Buried surface area (Å²) |
| SC | [0, 1] | Shape complementarity (standalone) |
| pLDDT | [0, 100] | Local prediction confidence |
| iPTM | [0, 1] | Interface PTM score |
| Clashes | [0, ∞] | Atomic steric clashes |
| Ramachandran | [0, 100]% | Favored region percentage |

### Shape Complementarity (SC)

Shape complementarity measures how well two protein surfaces fit together geometrically.
The current SC path is a standalone Python implementation inspired by Lawrence & Colman (1993)
J. Mol. Biol. 234: 946-950 and calibrated against the original CCP4 `sc` program.

Current implementation details:
- Uses CCP4-like parameter defaults (`DOT_DENSITY=15`, `TRIM=1.5`, `WEIGHT=0.5`, `INTERFACE=8`)
- Uses a copied CCP4-style `sc_radii.lib` atom-radius table
- Uses a Python SES/surface approximation rather than the original CCP4 Connolly surface engine
- Uses combined-molecule scoring for multi-chain partners in the active package path

Current limitations:
- **Binary two-chain interfaces are the most reliable use case today**
- **Three or more chains in one SC molecule are still limited**: the active package supports combined
  query/partner chain groups, but the current Python surface model is still an approximation and can
  deviate noticeably from CCP4, especially for multi-chain partner surfaces such as `A:(B+C)`
- Large structures can be slow and memory-intensive in the Python implementation
- The current Python SC implementation is **not yet a full replacement** for the original CCP4 Connolly
  surface construction on all datasets

Practical recommendation:
- For routine binary-interface scoring, the current Python SC path is usable
- For high-confidence validation, especially multi-chain SC or publication-grade comparison,
  cross-check against the original CCP4 `sc` program

### Preliminary-result caution

Current interface-analysis results should be treated **cautiously** for:
- multi-chain SC problems where one partner is composed of three or more chains
- membrane-protein systems

Reason:
- SC in the active package is still an approximation to the original CCP4 `sc` geometry engine
- multi-chain partner surfaces can behave differently from simple binary interfaces
- membrane-protein systems may stress both SC and legacy benchmarking assumptions

In practice:
- binary soluble interfaces are the most reliable current use case
- multimer and membrane-protein results should be treated as **preliminary** until further backend refinement and validation are completed

See `CURRENT_STATUS_SC_BSA.md` in the repository root for the current bottlenecks, benchmark notes,
and proposed next implementation step.

## Configuration

### Intercaat Configuration

Edit `intercaat_config.ini` to configure Voronoi tessellation:

```ini
[qvoronoi_path]
qvoronoi_bin = /path/to/qhull/bin/
executable_name = qvoronoi Fi
run_python_version = no  # Set to 'yes' to use scipy instead
```

### Chain Specification

The `-c` argument specifies which chains form the interface using a **hub-and-spoke model**:

#### Without separator (first chain is query):
- **First chain**: Query chain (main protein of interest)
- **Remaining chains**: Interaction chains (binding partners)

| Example | Description |
|---------|-------------|
| `-c A B` | Interface between chain A and chain B |
| `-c A B C` | Interface between chain A and chains B+C |
| `-c A B C D E` | Interface between chain A and chains B+C+D+E |

#### With `:` separator (multiple query chains):
Use `:` to explicitly separate query chains from partner chains.

| Example | Description |
|---------|-------------|
| `-c A : B` | Interface between chain A and chain B |
| `-c A : B C` | Interface between chain A and chains B+C |
| `-c A B : C D` | **A and B** each interface with **C and D** |
| `-c A B C : D E` | **A, B, C** each interface with **D and E** |

#### Example scenarios:

```bash
# Single query, multiple partners
# Complex like: A-B-C where A is the main protein
-c A B C

# Multiple queries, multiple partners
# Complex like: (A,B)-(C,D,E) where A,B are one subunit and C,D,E are another
-c A B : C D E

# Binary complex: A-B
-c A B
# or
-c A : B
```

**Note**: The interface metrics (contacts, clashes, BSA) are calculated between ALL query-partner pairs.

## MPI Parallelization

For large datasets, enable MPI to distribute work across cores:

```bash
# 4 processes
mpirun -np 4 int-iscore af3-ref -i /path/to/data -c A B --mpi

# Or with SLURM
srun -n 16 int-iscore md -i /path/to/trajectories -c A B --mpi
```

## Output

The tool generates a CSV file with columns:
- `filename`: Input file name
- `average_plddt`: Mean pLDDT score
- `ranking_score`: AlphaFold ranking metric
- `iptm/ptm`: Interface/full PTM scores
- `dope/soap`: Structure quality scores
- `dockq`: DockQ vs reference (af3-ref mode only)
- `sc`: Shape complementarity score
- `cpscore`: Contact preference score
- `bsa`: Buried surface area
- `interface_plddt`: Mean pLDDT at interface
- `ramachandran`: PROCHECK-style phi/psi analysis
- `interface_residues`: List of interface residue IDs

## Legacy utility wrappers

The package currently exposes two **legacy compatibility wrappers**:

- `int-iscore-frustration-legacy`
- `int-iscore-input-data-legacy`

These wrappers are intended to preserve access to historical repo-root workflows.
They are **not** yet fully package-native modes and should be treated as source-tree-dependent
compatibility entrypoints rather than stable release-grade APIs.

Example execution commands:

```bash
# Legacy frustratometer batch workflow
int-iscore-frustration-legacy

# Legacy benchmark / input-data generation workflow
int-iscore-input-data-legacy
```

These commands are expected to be run from a checked-out source tree that still contains
the historical supporting files and templates used by the original scripts.

## DockQ Quality Thresholds

| DockQ Range | Quality |
|-------------|---------|
| 0.00 - 0.23 | Incorrect |
| 0.23 - 0.49 | Acceptable |
| 0.49 - 0.80 | Medium |
| >= 0.80 | High |

## Ramachandran Quality Thresholds

- **> 95% in favored regions**: Good model
- **90-95%**: Acceptable
- **< 90%**: Poor geometry

## Citation

If you use int-iScore in your research, please cite:

```
Luo, Y. (2024). int-iScore: Integrated AI Predictions Interface Analysis Scores Toolkit.
GitHub Repository. https://github.com/lluoto/int-iScore
```

## License

MIT License

## Contact

**Author**: Yulin Luo  
**Email**: luoyl2022@shanghaitech.edu.cn  
**Institution**: ShanghaiTech University

## Acknowledgments

This toolkit integrates several excellent tools:
- [AlphaFold3](https://github.com/google-deepmind/alphafold3)
- [DockQ](https://github.com/bjornwallner/DockQ)
- [CPscore](https://github.com/wallnerlab/ProQDock)
- [Intercaat](https://github.com/eved1018/Intercaat)
- [MODELLER](https://salilab.org/modeller/)
- [freesasa](https://github.com/freesasa/freesasa-python)
