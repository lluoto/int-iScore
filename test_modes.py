#!/usr/bin/env python3
"""Test int-iScore modes on server."""
import sys
import os
import argparse

sys.path.insert(0, '/media/cuixi/data0/lluoto/int-iScore/src')
os.chdir('/media/cuixi/data0/lluoto/int-iScore')

def test_md():
    """Test md mode with split_frames (ABCD:E)."""
    from int_iscore.modes.md_snapshots import run_md_analysis
    args = argparse.Namespace(
        input="/media/cuixi/data0/lluoto/split_frames",
        output="/media/cuixi/data0/lluoto/int-iScore/output_md",
        name="md_test",
        chains=["A", "B", "C", "D", ":", "E"],
        cluster_method="spectral",
        n_clusters=5,
        cache_dir="cache",
        ccp4_path=None,
        mpi=False,
    )
    print("=== Testing MD mode (split_frames, ABCD:E) ===")
    run_md_analysis(args)

def test_af3_nonref():
    """Test af3-nonref mode with 5-HT2C (A:BC)."""
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
    run_nonref_analysis(args)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.makedirs("/media/cuixi/data0/lluoto/int-iScore/output_md", exist_ok=True)
    os.makedirs("/media/cuixi/data0/lluoto/int-iScore/output_5HT2C", exist_ok=True)
    
    if mode in ("md", "all"):
        test_md()
    if mode in ("af3-nonref", "all"):
        test_af3_nonref()
