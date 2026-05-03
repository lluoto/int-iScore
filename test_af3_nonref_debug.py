#!/usr/bin/env python3
"""Debug test for af3-nonref mode."""
import sys
import os
import argparse

sys.path.insert(0, '/media/cuixi/data0/lluoto/int-iScore/src')
os.chdir('/media/cuixi/data0/lluoto/int-iScore')

from int_iscore.modes.af3_nonref import run_nonref_analysis

args = argparse.Namespace(
    input="/media/cuixi/data0/lluoto/5-HT2C-2CD28",
    output="/media/cuixi/data0/lluoto/int-iScore/output_5HT2C",
    name="5HT2C_test",
    chains=["A", ":", "B", "C"],
    format="cif",
    temp_dir="int_iscore_temp",
    ccp4_path=None,
)

print("=== Testing AF3 non-ref mode (5-HT2C, A:BC) ===")
print("Input:", args.input)
print("Output:", args.output)
print("Chains:", args.chains)

# List input files
cif_files = [f for f in os.listdir(args.input) if f.endswith('.cif')]
print(f"Found {len(cif_files)} CIF files: {cif_files[:3]}...")

run_nonref_analysis(args)
