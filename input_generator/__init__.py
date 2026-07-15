"""
input_generator - General-purpose protein structure prediction input generator.

AF3/Chai/Boltz input generation with intelligent sequence trimming and MSA support.

Features:
  - Sequence trimming via 3 methods (PDB CA, UniProt TM, Kyte-Doolittle)
  - AF3 (alphafoldserver dialect) input JSON generation
  - Boltz-1 input generation
  - Chai-1 input generation
  - ColabFold MMseqs2 MSA computation
  - Pair-wise batch processing (homo/heterodimers)
"""

from .trimming import trim_to_7tm, load_pdb_ca_ranges, load_uniprot_tm_data
from .trimming import predict_tm_regions, get_7tm_boundaries
from .msa import compute_msa, parse_a3m_sequences
from .generators import (
    make_af3_json,
    make_boltz_input,
    make_chai_input,
    generate_all_formats,
)

__all__ = [
    # Trimming
    'trim_to_7tm', 'load_pdb_ca_ranges', 'load_uniprot_tm_data',
    'predict_tm_regions', 'get_7tm_boundaries',
    # MSA
    'compute_msa', 'parse_a3m_sequences',
    # Generators
    'make_af3_json', 'make_boltz_input', 'make_chai_input',
    'generate_all_formats',
]
