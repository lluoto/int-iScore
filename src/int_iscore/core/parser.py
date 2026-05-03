"""
Command-line argument parsing and configuration utilities.
"""

import argparse
import os
from pathlib import Path


def parse_config(config_file="intercaat_config.ini"):
    """
    Parse intercaat configuration file.
    
    Args:
        config_file: Path to intercaat config file
    
    Returns:
        Dictionary with configuration values
    """
    from configparser import ConfigParser
    
    config = ConfigParser()
    config.read(config_file)
    
    return {
        "qvoronoi_bin": config.get("qvoronoi_path", "qvoronoi_bin", fallback=""),
        "executable_name": config.get("qvoronoi_path", "executable_name", fallback="qvoronoi Fi"),
        "run_python_version": config.get("qvoronoi_path", "run_python_version", fallback="no"),
    }


def create_base_parser():
    """Create base argument parser with common arguments."""
    parser = argparse.ArgumentParser(
        description="int-iScore: Protein complex interface analysis toolkit",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input path (directory or file)",
        type=str,
    )
    parser.add_argument(
        "-o", "--output",
        default=".",
        help="Output directory for results",
        type=str,
    )
    parser.add_argument(
        "-n", "--name",
        default="result",
        help="Output file prefix",
        type=str,
    )
    parser.add_argument(
        "-c", "--chains",
        default=["A", "B"],
        nargs="+",
        help="Chain identifiers for interface analysis. "
             "Use ':' to separate query chains from partner chains. "
             "Example: 'A B : C D E' means chains A,B query, chains C,D,E are partners. "
             "Example: 'A : B C D' means chain A queries, chains B,C,D are partners.",
        type=str,
    )
    parser.add_argument(
        "--bsa-method",
        choices=["freesasa", "vmd"],
        default="freesasa",
        help="Method for buried surface area calculation",
    )
    parser.add_argument(
        "--clash-cutoff",
        default=0.4,
        type=float,
        help="Clash detection threshold multiplier",
    )
    
    return parser


def parse_chain_specification(chains_arg):
    """
    Parse chain specification with separator.
    
    Supports two formats:
    1. Without separator: ['A', 'B', 'C'] -> queries=['A'], partners=['B', 'C']
    2. With separator: ['A', ':', 'B', 'C'] -> queries=['A'], partners=['B', 'C']
    3. With separator: ['A', 'B', ':', 'C', 'D'] -> queries=['A', 'B'], partners=['C', 'D']
    
    Args:
        chains_arg: List of chain strings from argparse
        
    Returns:
        Tuple of (query_chains, partner_chains)
    """
    if isinstance(chains_arg, str):
        chains_arg = chains_arg.split()
    
    if ":" in chains_arg:
        sep_idx = chains_arg.index(":")
        queries = [c for c in chains_arg[:sep_idx] if c != ":"]
        partners = [c for c in chains_arg[sep_idx+1:] if c != ":"]
        return queries, partners
    else:
        if len(chains_arg) >= 2:
            return [chains_arg[0]], chains_arg[1:]
        return chains_arg, []


def parse_af3_args():
    """Parse arguments for AlphaFold3 reference mode."""
    parser = create_base_parser()
    
    parser.description = "AlphaFold3 with reference - Calculate DockQ and interface metrics"
    parser.add_argument(
        "-r", "--reference",
        help="Reference PDB file (optional, will auto-detect if not provided)",
        type=str,
    )
    parser.add_argument(
        "--reference-dir",
        default="reference",
        help="Directory containing reference PDB files",
        type=str,
    )
    parser.add_argument(
        "--server-output",
        action="store_true",
        help="Use full_data JSON instead of confidences for pLDDT",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Enable MPI parallel processing",
    )
    parser.add_argument(
        "--no-clean-pdb",
        action="store_true",
        help="Skip automatic PDB cleaning",
    )
    
    return parser.parse_args()


def parse_md_args():
    """Parse arguments for molecular dynamics snapshots mode."""
    parser = create_base_parser()
    
    parser.description = "Molecular dynamics snapshots - Cluster and analyze trajectory frames"
    parser.add_argument(
        "--cluster-method",
        choices=["spectral", "hierarchical", "dbscan"],
        default="spectral",
        help="Clustering method for representative selection",
    )
    parser.add_argument(
        "--n-clusters",
        default=20,
        type=int,
        help="Number of clusters for spectral clustering",
    )
    parser.add_argument(
        "--rmsd-cutoff",
        default=0.5,
        type=float,
        help="RMSD cutoff for DBSCAN clustering",
    )
    parser.add_argument(
        "--min-samples",
        default=5,
        type=int,
        help="Minimum samples for DBSCAN",
    )
    parser.add_argument(
        "--cache-dir",
        default="cache",
        help="Directory for cached RMSD matrices",
        type=str,
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Enable MPI parallel processing",
    )
    parser.add_argument(
        "--chimera-path",
        help="Path to Chimera executable or installation for EC calculation",
        type=str,
    )
    parser.add_argument(
        "--vmd-path",
        help="Path to VMD executable for EC calculation",
        type=str,
    )
    parser.add_argument(
        "--pdb2pqr-path",
        help="Path to pdb2pqr executable for EC calculation",
        type=str,
    )
    parser.add_argument(
        "--apbs-path",
        help="Path to APBS executable for EC calculation",
        type=str,
    )
    
    return parser.parse_args()


def parse_nonref_args():
    """Parse arguments for non-reference AlphaFold3 mode."""
    parser = create_base_parser()
    
    parser.description = "AlphaFold3 without reference - Interface metrics only"
    parser.add_argument(
        "--format",
        choices=["auto", "cif", "pdb"],
        default="auto",
        help="Input file format; auto selects any supported structure files present",
    )
    parser.add_argument(
        "--temp-dir",
        default="int_iscore_temp",
        help="Directory for temporary converted files",
        type=str,
    )
    parser.add_argument(
        "--chimera-path",
        help="Path to Chimera executable or installation for EC calculation",
        type=str,
    )
    parser.add_argument(
        "--vmd-path",
        help="Path to VMD executable for EC calculation",
        type=str,
    )
    parser.add_argument(
        "--pdb2pqr-path",
        help="Path to pdb2pqr executable for EC calculation",
        type=str,
    )
    parser.add_argument(
        "--apbs-path",
        help="Path to APBS executable for EC calculation",
        type=str,
    )
    
    return parser.parse_args()


def parse_single_args():
    """Parse arguments for single structure analysis mode."""
    parser = create_base_parser()
    
    parser.description = "Single structure analysis - Quick interface metrics"
    parser.add_argument(
        "input",
        help="Input PDB/CIF file",
        type=str,
    )
    parser.add_argument(
        "-o", "--output",
        default="single_result.csv",
        help="Output CSV file",
        type=str,
    )
    parser.add_argument(
        "-c", "--chains",
        default=["A", "B"],
        nargs="+",
        help="Chain identifiers",
        type=str,
    )
    
    return parser.parse_args()
