#!/usr/bin/env python
"""
int-iScore: Integrated Interface Scores Toolkit

Command-line interface for protein complex interface analysis.
"""

import argparse
import sys
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="int-iScore: Integrated Interface Scores Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # AlphaFold3 with reference (DockQ calculation)
  int-iscore af3-ref -i /path/to/af3/output -c A B -n my_analysis

  # Molecular dynamics snapshots (with clustering)
  int-iscore md -i /path/to/md/trajectories -c A B --cluster-method spectral

  # AlphaFold3 without reference
  int-iscore af3-nonref -i /path/to/af3/output -c A B -n my_analysis

Modes:
  af3-ref     AlphaFold3 models with reference structure (calculates DockQ)
  md          Molecular dynamics trajectory snapshots (with RMSD clustering)
  af3-nonref  AlphaFold3 models without reference structure
        """
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="int-iScore 1.0.0"
    )
    
    subparsers = parser.add_subparsers(dest="mode", help="Analysis mode")
    
    parser_af3_ref = subparsers.add_parser(
        "af3-ref",
        help="AlphaFold3 with reference structure"
    )
    parser_af3_ref.add_argument("-i", "--input", required=True, help="Input directory")
    parser_af3_ref.add_argument("-c", "--chains", nargs="+", default=["A", "B"],
                                help="Chain identifiers: first is query chain, rest are partners. "
                                     "Example: 'A B C' means A interfaces with B and C. "
                                     "Use ':' separator: 'A B : C D' means A,B query with C,D partners")
    parser_af3_ref.add_argument("-n", "--name", default="result", help="Output prefix")
    parser_af3_ref.add_argument("-o", "--output", default=".", help="Output directory")
    parser_af3_ref.add_argument("-r", "--reference", help="Reference PDB file")
    parser_af3_ref.add_argument("--reference-dir", default="reference",
                                help="Reference directory")
    parser_af3_ref.add_argument("--server-output", action="store_true",
                                help="Use full_data JSON")
    parser_af3_ref.add_argument("--ccp4-path", help="Path to CCP4 installation for SC calculation")
    parser_af3_ref.add_argument("--chimera-path", help="Path to Chimera executable or installation for EC calculation")
    parser_af3_ref.add_argument("--vmd-path", help="Path to VMD executable for EC calculation")
    parser_af3_ref.add_argument("--pdb2pqr-path", help="Path to pdb2pqr executable for EC calculation")
    parser_af3_ref.add_argument("--apbs-path", help="Path to APBS executable for EC calculation")
    parser_af3_ref.add_argument("--mpi", action="store_true", help="Enable MPI")
    parser_af3_ref.add_argument("--no-clean-pdb", action="store_true", help="Skip automatic PDB cleaning")
    
    parser_md = subparsers.add_parser(
        "md",
        help="Molecular dynamics snapshots"
    )
    parser_md.add_argument("-i", "--input", required=True, help="Input directory")
    parser_md.add_argument("-c", "--chains", nargs="+", default=["A", "B"],
                           help="Chain identifiers: first is query chain, rest are partners. "
                                "Example: 'A B C' means A interfaces with B and C. "
                                "Use ':' separator: 'A B : C D' means A,B query with C,D partners")
    parser_md.add_argument("-n", "--name", default="result", help="Output prefix")
    parser_md.add_argument("-o", "--output", default=".", help="Output directory")
    parser_md.add_argument("--cluster-method", choices=["spectral", "hierarchical", "dbscan"],
                          default="spectral", help="Clustering method")
    parser_md.add_argument("--n-clusters", type=int, default=20, help="Number of clusters")
    parser_md.add_argument("--cache-dir", default="cache", help="Cache directory")
    parser_md.add_argument("--ccp4-path", help="Path to CCP4 installation for SC calculation")
    parser_md.add_argument("--chimera-path", help="Path to Chimera executable or installation for EC calculation")
    parser_md.add_argument("--vmd-path", help="Path to VMD executable for EC calculation")
    parser_md.add_argument("--pdb2pqr-path", help="Path to pdb2pqr executable for EC calculation")
    parser_md.add_argument("--apbs-path", help="Path to APBS executable for EC calculation")
    parser_md.add_argument("--mpi", action="store_true", help="Enable MPI")
    
    parser_af3_nonref = subparsers.add_parser(
        "af3-nonref",
        help="AlphaFold3 without reference"
    )
    parser_af3_nonref.add_argument("-i", "--input", required=True, help="Input directory")
    parser_af3_nonref.add_argument("-c", "--chains", nargs="+", default=["A", "B"],
                                    help="Chain identifiers: first is query chain, rest are partners. "
                                         "Example: 'A B C' means A interfaces with B and C. "
                                         "Use ':' separator: 'A B : C D' means A,B query with C,D partners")
    parser_af3_nonref.add_argument("-n", "--name", default="result", help="Output prefix")
    parser_af3_nonref.add_argument("-o", "--output", default=".", help="Output directory")
    parser_af3_nonref.add_argument("--format", choices=["auto", "cif", "pdb"], default="auto",
                                   help="Input format; auto selects any supported structure files present")
    parser_af3_nonref.add_argument("--temp-dir", default="int_iscore_temp",
                                   help="Temporary directory")
    parser_af3_nonref.add_argument("--ccp4-path", help="Path to CCP4 installation for SC calculation")
    parser_af3_nonref.add_argument("--chimera-path", help="Path to Chimera executable or installation for EC calculation")
    parser_af3_nonref.add_argument("--vmd-path", help="Path to VMD executable for EC calculation")
    parser_af3_nonref.add_argument("--pdb2pqr-path", help="Path to pdb2pqr executable for EC calculation")
    parser_af3_nonref.add_argument("--apbs-path", help="Path to APBS executable for EC calculation")
    
    args = parser.parse_args()
    
    if args.mode is None:
        parser.print_help()
        sys.exit(1)
    
    # Save original working directory
    original_cwd = os.getcwd()
    SCRIPT_DIR = Path(__file__).parent
    os.chdir(SCRIPT_DIR)
    
    # Pass original working directory to mode functions
    args._original_cwd = original_cwd
    
    if args.mode == "af3-ref":
        from int_iscore.modes.af3_reference import run_af3_analysis
        run_af3_analysis(args)
    elif args.mode == "md":
        from int_iscore.modes.md_snapshots import run_md_analysis
        run_md_analysis(args)
    elif args.mode == "af3-nonref":
        from int_iscore.modes.af3_nonref import run_nonref_analysis
        run_nonref_analysis(args)


if __name__ == "__main__":
    main()
