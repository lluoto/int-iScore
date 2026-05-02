"""
Shape Complementarity (SC) calculation module.

Uses standalone Python implementation based on Lawrence & Colman (1993).
This does not require CCP4 and can handle large protein structures.

Reference:
- Lawrence & Colman (1993) Shape complementarity at protein-protein interfaces
  J. Mol. Biol. 234: 946-950.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional

from .sc_calculator import (
    read_pdb_atoms,
    calculate_sc_from_pdb,
    calculate_sc_combined_from_pdb,
    calculate_sc_batch as calc_sc_batch,
    calculate_interface_atoms,
    PROBE_RADIUS,
    DOT_DENSITY,
    WEIGHT,
)


def calculate_sc(
    pdb_file: str,
    query_chains: List[str],
    partner_chains: List[str],
    probe_radius: float = PROBE_RADIUS,
    dot_density: int = DOT_DENSITY,
    weight: float = WEIGHT,
    interface_distance: float = 8.0,
) -> Dict:
    """
    Calculate SC scores for a protein complex.
    
    Args:
        pdb_file: Path to PDB file
        query_chains: List of query chain IDs
        partner_chains: List of partner chain IDs
        probe_radius: Radius of probe sphere for surface calculation
        dot_density: Density of surface points per A^2
        weight: Weighting factor for S(A->B) calculation
        interface_distance: Distance to consider atoms at interface
    
    Returns:
        Dictionary with SC scores per chain pair and statistics
    """
    combined_score = calculate_sc_combined_from_pdb(
        pdb_file,
        query_chains,
        partner_chains,
        probe_radius=probe_radius,
        dot_density=dot_density,
        weight=weight,
        interface_distance=interface_distance,
    )
    sc_scores = {("combined_query", "combined_partner"): combined_score}

    if len(query_chains) == 1 and len(partner_chains) == 1:
        sc_scores = calc_sc_batch(
            pdb_file,
            query_chains,
            partner_chains,
            probe_radius=probe_radius,
            dot_density=dot_density,
            weight=weight,
            interface_distance=interface_distance,
        )
    
    results = {
        'sc_scores': sc_scores,
        'sc_mean': 0.0,
        'sc_min': 0.0,
        'sc_max': 0.0,
        'n_pairs': len(sc_scores),
    }
    
    if sc_scores:
        score_values = list(sc_scores.values())
        results['sc_mean'] = float(np.mean(score_values))
        results['sc_min'] = float(np.min(score_values))
        results['sc_max'] = float(np.max(score_values))
    
    return results


def calculate_sc_pair(
    pdb_file: str,
    chain1: str,
    chain2: str,
    **kwargs
) -> float:
    """
    Calculate SC for a single chain pair.
    
    Args:
        pdb_file: Path to PDB file
        chain1: First chain ID
        chain2: Second chain ID
        **kwargs: Additional arguments for SC calculation
    
    Returns:
        SC score (0 to 1, higher is better)
    """
    return calculate_sc_from_pdb(pdb_file, chain1, chain2, **kwargs)


def get_interface_atoms(
    pdb_file: str,
    chain1: str,
    chain2: str,
    interface_distance: float = 8.0
) -> Tuple[List[Dict], List[Dict]]:
    """
    Get interface atoms between two chains.
    
    Args:
        pdb_file: Path to PDB file
        chain1: First chain ID
        chain2: Second chain ID
        interface_distance: Distance threshold for interface atoms
    
    Returns:
        Tuple of (interface_atoms_chain1, interface_atoms_chain2)
    """
    atoms1 = read_pdb_atoms(pdb_file, [chain1])
    atoms2 = read_pdb_atoms(pdb_file, [chain2])
    return calculate_interface_atoms(atoms1, atoms2, interface_distance)


def check_sc_available() -> bool:
    """
    Check if SC calculation is available (always True for standalone).
    
    Returns:
        Always True for standalone implementation
    """
    return True
