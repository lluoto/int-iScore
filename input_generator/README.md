"""
input_generator - General-purpose protein structure prediction input generator.

AF3/Chai/Boltz input generation with intelligent sequence trimming and MSA support.

## Features

- **Sequence trimming** with 3-tier priority:
  1. UniProt TRANSMEM annotations (curated 7TM boundaries)
  2. PDB CA atom coordinates (actual resolved residues)
  3. Kyte-Doolittle hydropathy sliding window (algorithmic fallback)

- **MSA computation** via ColabFold MMseqs2 server
- **Multi-format output**: AF3 (alphafoldserver dialect), Boltz-1, Chai-1
- **Homodimer/heterodimer/monomer** support
- **Batch processing** with caching and resume

## Quick Start

```python
from input_generator import trim_sequence, generate_all_formats

# Trim sequence
trimmed, start, end, method = trim_sequence(
    seq='MVLSPADK...',
    uniprot_id='P07550',
    pdb_ca_path='pdb_ca_ranges.json',
    uniprot_tm_path='uniprot_tm_boundaries.json',
)
print('Trimmed: {} aa, method: {}'.format(len(trimmed), method))

# Generate all model inputs
files = generate_all_formats(
    name='my_protein',
    seq1=trimmed,
    output_dir='./output',
    seeds=20,
)
print('Generated:', files)
```

## Module Structure

```
input_generator/
  __init__.py     # Package exports
  trimming.py     # Sequence trimming (PDB/UniProt/Kyte-Doolittle)
  msa.py          # MSA computation via ColabFold
  generators.py   # AF3/Boltz/Chai input JSON generation
```

## Dependencies

- Python 3.8+
- BioPython (optional, for PDB parsing in trimming)
- requests (for ColabFold MSA server)
- mmseqs2 (from ColabFold, for MSA computation)

## Integration with int-iScore

Clone this module into the int-iScore repo:
```bash
cp -r input_generator /path/to/int-iscore/
```

Then use alongside the scoring pipeline:
```python
from input_generator import trim_sequence, generate_all_formats

# Step 1: Trim & generate inputs
generate_all_formats(...)

# Step 2: Predict with AF3/Boltz/Chai

# Step 3: Score with int-iScore
# (use the batch_results/ scoring pipeline)
```
