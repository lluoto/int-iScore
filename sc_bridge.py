#!/usr/bin/env python3
"""Bridge script: called by sc_hts to get CCP4-accurate SC score."""
import sys
import os

if len(sys.argv) != 4:
    print("Usage: sc_bridge.py <pdb_path> <chain1> <chain2>")
    sys.exit(1)

pdb_path = sys.argv[1]
chain1 = sys.argv[2]
chain2 = sys.argv[3]

# Ensure we are in the project root directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, "src")
from int_iscore.utils.sc_calculator import calculate_sc_from_pdb

try:
    sc = calculate_sc_from_pdb(pdb_path, chain1, chain2)
    print(sc)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
