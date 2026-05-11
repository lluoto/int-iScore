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

## SC Optimization (2026-05-09)

The Shape Complementarity (SC) calculator has been significantly improved:

- **Surface Generation**: Integrated SCASA connolly.py (CCP4 mds Fortran-to-Python
  translation, GPL-3.0) for full Connolly molecular surface generation (convex +
  toroidal + concave patches). Replaces the previous SAS (solvent-accessible surface)
  approximation.
- **Scoring Formula**: Fixed to match CCP4's negate convention: `-(n1.n2) * exp(-w * d^2)`.
- **NN Search Target**: Changed from full surface to trimmed surface (PiA->PiB),
  matching the original Lawrence & Colman (1993) algorithm.

### Verification against CCP4 Reference

| PDB | Chains | SC (int-iScore) | SC (CCP4) | Delta |
|-----|--------|----------------|-----------|-------|
| 6A6I | A-B | 0.607 | 0.616 | -1.5% |
| 5HT2C | A-B | 0.453 | 0.448 | +1.1% |

Both results are within the expected variation of the SCASA reference implementation
(delta <0.02 from CCP4).

### SCASA License

The Connolly molecular surface generation (`connolly.py`) is derived from the
[SCASA](https://github.com/WhalleyT/SCASA) project and is distributed under the
**GNU General Public License v3.0 (GPL-3.0)**. See `SCASA_LICENSE` for details.
```

### C++ High-Throughput SC Engine (sc_hts)

A C++17 two-tier Shape Complementarity engine for screening 10,000+ protein complexes.

**Architecture:**
- **Mode 0 (FastSAS, EXPERIMENTAL)**: Two scoring methods (Voxel density + Feature regression). ~200 ms/complex. **Under active development — does NOT yet generalize across PDBs.** See FASTSAS_ANALYSIS_PROMPT.md for current status.
  Coarse filter for ranking large batches. Gaussian Overlap scoring (normal-vector-dominant): weight=0.05, distance cutoff=4.0 Angstrom.
- **Mode 1 (Accurate, DEFAULT)**: Calls Python SCASA bridge via popen for CCP4-accurate SC (±2%)
  (within +/-2%). ~9-50 s/complex depending on complex size.
- **Mode 2 (Batch)**: Reads CSV job list, processes sequentially, outputs combined CSV.

**Build:**
```bash
g++ -std=c++17 -O3 -fopenmp -I. sc_hts.cpp -o sc_hts
```

**Usage:**
```bash
# Single PDB, FastSAS (coarse)
./sc_hts complex.pdb A B 0 pdb_id

# Single PDB, Accurate (Python SCASA)
./sc_hts complex.pdb A B 1 pdb_id

# Batch mode from CSV (pdb_path,chain1,chain2,mode[,pdb_id])
./sc_hts jobs.csv ? ? 2
```

**Benchmark (6A6I A-B):**

| Mode      | SC      | Time     | Dots   | Accuracy       |
|-----------|---------|----------|--------|----------------|
| FastSAS (experimental) | 0.58    | 200 ms   | 400K   | Under development (does not generalize) |
| Accurate  | 0.585   | 9.5 s    | -      | CCP4 +/-2%   |
| CCP4 ref  | 0.616   | -        | -      | Gold standard  |

**FastSAS Development Status:** FastScoring is under active development. Current approaches:
- **Voxel density** (v2): poor cross-PDB discrimination (r = -0.34)
- **4-feature linear regression** (v3): overfits training set (calibration r=+0.75, cross-validation r=-0.15)
- **GBDT dual-track model** (v4): planned, see FASTSAS_ANALYSIS_PROMPT.md

**For production use, Accurate mode is the default and only reliable option.**

**Dependencies:** g++ 9+, OpenMP, nanoflann.hpp (bundled), Python 3 with int_iscore package


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
- `frustration_score`: Configurational frustration (Frustratometer2) sum over interface
- `bsa`: Buried surface area
- `interface_plddt`: Mean pLDDT at interface
- `ramachandran`: PROCHECK-style phi/psi analysis
- `interface_residues`: List of interface residue IDs

## Frustration Scoring (Frustratometer2)

The toolkit integrates [Frustratometer2](https://github.com/proteinphysiologylab/frustratometer2)
for configurational frustration analysis of protein-protein interfaces. Frustration scores
quantify the energetic local frustration at residue-residue contacts based on a
AWSEM/LAMMPS statistical mechanics pipeline.

### Integration

- **Active pipeline**: `md_snapshots.py` computes `frustration_score` for each
  trajectory snapshot by calling `core.metrics.compute_frustration_score()`
- **Legacy wrapper**: `int-iscore-frustration-legacy` preserves the original
  frustratometer2 batch workflow for standalone use

### Method

For each interface residue pair, Frustratometer2 computes a Frustration Index (FrstIndex)
comparing the native (or predicted) energy to a decoy ensemble. The `frustration_score`
column reports the sum of all inter-chain FrstIndex values at the interface.

### Dependencies

- Perl 5+ (`/usr/bin/perl`)
- LAMMPS serial binary (`lmp_serial_12`) included in `frustratometer2-master/Scripts/`
- AWSEM data files (included)

### Usage

Frustration scoring in `md-snapshots` mode is automatic when the `frustratometer2-master/`
directory is present in the repository root. No additional configuration is needed.

## Legacy utility wrappers

The package also exposes two legacy compatibility wrappers for historical workflows:

- `int-iscore-frustration-legacy`
- `int-iscore-input-data-legacy`

These wrappers are **not** package-native modes and depend on the full source tree.
They are maintained for backward compatibility only.

Example:

```bash
# Legacy frustratometer batch workflow
int-iscore-frustration-legacy

# Legacy benchmark / input-data generation workflow
int-iscore-input-data-legacy
```

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
