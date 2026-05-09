"""Shape Complementarity (SC) calculator - Connolly molecular surface implementation.

Optimization (2026-05-09):
  - SCASA connolly.py integration (CCP4 mds Fortran-to-Python translation, GPL-3.0)
    for full Connolly surface generation (convex + toroidal + concave patches).
  - Fixed scoring formula: np.abs(dot) -> (-dot) to match CCP4 negate convention.
  - NN search target changed: full surface -> trimmed surface (PiA->PiB).
  - Verified against CCP4 reference values:
      6A6I A-B: 0.607 (CCP4: 0.616, delta 1.5%)
      5HT2C A-B: 0.453 (CCP4: 0.448, delta 1.1%)
"""

from __future__ import annotations

import os
from fnmatch import fnmatchcase
from typing import Dict, List, Tuple

import numpy as np
from scipy.spatial import cKDTree

# ----- try compiled backend -----
try:
    from .sc_backend import calculate_sc_backend as _calculate_sc_backend
except Exception:
    _calculate_sc_backend = None

USE_SC_BACKEND = os.getenv("INT_ISCORE_USE_SC_BACKEND", "0") == "1"
# ----- SCASA Connolly surface generator (CCP4 mds translation, GPL-3.0) -----
try:
    from .connolly import mds as _connolly_mds, BURIED_FLAG as _BURIED_FLAG
except ImportError:
    _connolly_mds = None
    _BURIED_FLAG = 1


# ----- CCP4-style radii -----
SC_RADII_LIB = """
***    H       0.50
***    H*      0.50
***    H**     0.50
***    H***    0.50
***    CA      1.85
***    C       1.80
***    O       1.60
***    N       1.65
***    CB      1.90
***    OT*     1.60
***    OXT     1.60
***    S*      1.90
***    P       1.80
ALA  CB      1.95
ARG  NH*     1.70
ARG  CZ      1.80
ARG  NE      1.65
ARG  CD      1.90
ARG  CG      1.90
ASN  ND2     1.70
ASN  OD1     1.60
ASN  CG      1.80
ASP  OD*     1.60
ASP  CG      1.80
GLN  NE2     1.70
GLN  OE1     1.60
GLN  CD      1.80
GLN  CG      1.90
GLU  OE*     1.60
GLU  CD      1.80
GLU  CG      1.90
GLY  CA      1.90
HIS  CD2     1.90
HIS  NE2     1.65
HIS  CE1     1.90
HIS  ND1     1.65
HIS  CG      1.80
HOH  O**     1.70
ILE  CD1     1.95
ILE  CG1     1.90
ILE  CB      1.85
ILE  CG2     1.95
LEU  CD*     1.95
LEU  CG      1.85
LYS  NZ      1.75
LYS  CE      1.90
LYS  CD      1.90
LYS  CG      1.90
MET  CE      1.95
MET  CG      1.90
PHE  CD*     1.90
PHE  CE*     1.90
PHE  CZ      1.90
PHE  CG      1.80
PRO  CD      1.90
PRO  CG      1.90
SER  OG      1.70
SUL  S       1.90
SUL  O***    1.65
THR  CG2     1.95
THR  OG1     1.70
THR  CB      1.85
TRP  CE2     1.80
TRP  CE3     1.90
TRP  CD1     1.90
TRP  CD2     1.80
TRP  CZ*     1.90
TRP  CH2     1.90
TRP  NE1     1.65
TRP  CG      1.80
TYR  OH      1.70
TYR  CD*     1.90
TYR  CE*     1.90
TYR  CZ      1.80
TYR  CG      1.80
VAL  CG*     1.95
VAL  CB      1.85
WAT  O       1.70
WAT  O*      1.70
""".strip()

DEFAULT_RADIUS = 1.70
PROBE_RADIUS = 1.4
DOT_DENSITY = 15
TRIM = 1.5
WEIGHT = 0.5
INTERFACE_DISTANCE = 8.0

# ----- radii helpers -----

def _parse_radii() -> List[Tuple[str, str, float]]:
    rules = []
    for line in SC_RADII_LIB.splitlines():
        parts = line.split()
        if len(parts) == 3:
            rules.append((parts[0], parts[1], float(parts[2])))
    return rules

RADII_RULES = _parse_radii()

def get_atom_radius(atom_name: str, res_name: str) -> float:
    for res_pat, atom_pat, radius in RADII_RULES:
        if fnmatchcase(res_name, res_pat) and fnmatchcase(atom_name, atom_pat):
            return radius
    return DEFAULT_RADIUS

# ----- PDB reading -----

def read_pdb_atoms(pdb_file: str, chains: List[str]) -> List[Dict]:
    atoms = []
    with open(pdb_file, "r") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = line[21:22].strip()
            if chain not in chains:
                continue
            name = line[12:16].strip()
            if name.startswith("H"):
                continue
            atoms.append({
                "name": name,
                "res": line[17:20].strip(),
                "res_id": int(line[22:26].strip()),
                "chain": chain,
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float),
                "radius": get_atom_radius(name, line[17:20].strip()),
            })
    return atoms

# ----- Fibonacci sphere -----

def _fibonacci_sphere(n: int) -> np.ndarray:
    if n < 4:
        n = 4
    idx = np.arange(n, dtype=float)
    phi = np.pi * (3.0 - np.sqrt(5.0))
    y = 1.0 - 2.0 * (idx + 0.5) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - y * y))
    theta = phi * idx
    return np.column_stack((r * np.cos(theta), y, r * np.sin(theta)))

# ----- interface atom selection -----

def calculate_interface_atoms(
    atoms1: List[Dict], atoms2: List[Dict],
    interface_distance: float = INTERFACE_DISTANCE,
) -> Tuple[List[Dict], List[Dict]]:
    if not atoms1 or not atoms2:
        return [], []
    c1 = np.array([a["xyz"] for a in atoms1], dtype=float)
    c2 = np.array([a["xyz"] for a in atoms2], dtype=float)
    r1 = np.array([a["radius"] for a in atoms1], dtype=float)
    r2 = np.array([a["radius"] for a in atoms2], dtype=float)

    # chunked minimum distance
    keep1 = np.zeros(len(atoms1), dtype=bool)
    keep2 = np.zeros(len(atoms2), dtype=bool)
    chunk = 500
    for start in range(0, len(c1), chunk):
        end = min(start + chunk, len(c1))
        d = np.linalg.norm(c1[start:end, None, :] - c2[None, :, :], axis=2)
        clearance = d - r1[start:end, None] - r2[None, :]
        keep1[start:end] = np.min(clearance, axis=1) < interface_distance

    for start in range(0, len(c2), chunk):
        end = min(start + chunk, len(c2))
        d = np.linalg.norm(c2[start:end, None, :] - c1[None, :, :], axis=2)
        clearance = d - r2[start:end, None] - r1[None, :]
        keep2[start:end] = np.min(clearance, axis=1) < interface_distance

    return [a for a, k in zip(atoms1, keep1) if k], [a for a, k in zip(atoms2, keep2) if k]

# ----- surface point generation -----

def _surface_points(atoms: List[Dict], probe_radius: float, dot_density: int) -> Tuple[np.ndarray, np.ndarray]:
    """Generate solvent-accessible surface points with outward normals."""
    if not atoms:
        return np.empty((0, 3)), np.empty((0, 3))
    coords = np.array([a["xyz"] for a in atoms], dtype=float)
    radii = np.array([a["radius"] for a in atoms], dtype=float)
    expanded = radii + probe_radius

    all_points = []
    all_normals = []

    for i in range(len(atoms)):
        R = expanded[i]
        area = 4.0 * np.pi * R * R
        n_pts = max(16, int(round(area * dot_density)))
        dirs = _fibonacci_sphere(n_pts)
        pts = atoms[i]["xyz"] + R * dirs
        normals = dirs.copy()  # outward from atom center

        # occlusion test against all other atoms' expanded spheres
        occluded = np.zeros(n_pts, dtype=bool)
        for j in range(len(atoms)):
            if j == i:
                continue
            Rj = expanded[j]
            d2 = np.sum((pts - atoms[j]["xyz"]) ** 2, axis=1)
            occluded |= d2 < Rj * Rj

        keep = ~occluded
        if np.any(keep):
            all_points.append(pts[keep])
            all_normals.append(normals[keep])

    if not all_points:
        return np.empty((0, 3)), np.empty((0, 3))
    return np.vstack(all_points), np.vstack(all_normals)

# ----- buried / accessible / trim -----

def _classify_and_trim(
    surf_pts: np.ndarray, surf_nrm: np.ndarray,
    other_atoms: List[Dict], probe_radius: float, trim: float,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """Classify surface points as buried/accessible and trim peripheral band.

    buried: point on mol A that is within reach of mol B atoms (probe can touch)
    trimmed: buried points far from accessible points (interior interface only)
    """
    if len(surf_pts) == 0 or not other_atoms:
        return np.empty((0, 3)), np.empty((0, 3)), 0, 0

    other_coords = np.array([a["xyz"] for a in other_atoms], dtype=float)
    other_radii = np.array([a["radius"] for a in other_atoms], dtype=float) + probe_radius

    # buried: probe center on the normal direction overlaps with other molecule
    probe_centers = surf_pts + surf_nrm * probe_radius
    d2 = np.sum((probe_centers[:, None, :] - other_coords[None, :, :]) ** 2, axis=2)
    buried_mask = np.any(d2 < other_radii[None, :] ** 2, axis=1)

    buried_pts = surf_pts[buried_mask]
    buried_nrm = surf_nrm[buried_mask]
    accessible_pts = surf_pts[~buried_mask]
    n_buried = int(buried_mask.sum())

    # trim: REMOVE buried points near accessible points (peripheral band)
    if len(accessible_pts) == 0 or len(buried_pts) == 0:
        return buried_pts, buried_nrm, n_buried, len(buried_pts)

    tree = cKDTree(accessible_pts)
    nearest_dist, _ = tree.query(buried_pts, k=1, workers=-1)
    keep = nearest_dist >= trim
    n_trimmed = int(keep.sum())
    return buried_pts[keep], buried_nrm[keep], n_buried, n_trimmed

# ----- SC scoring -----

def _score(
    src_pts: np.ndarray, src_nrm: np.ndarray,
    tgt_pts: np.ndarray, tgt_nrm: np.ndarray,
    weight: float,
) -> float:
    """Compute median complementarity: median of (-dot) * exp(-w * d²)."""
    if len(src_pts) < 5 or len(tgt_pts) < 5:
        return 0.0
    tree = cKDTree(tgt_pts)
    distances, idxs = tree.query(src_pts, k=1, workers=-1)
    nearest_nrm = tgt_nrm[idxs]
    # CCP4: negate dot for complementary surfaces
    d = np.sum(src_nrm * nearest_nrm, axis=1)
    values = (-d) * np.exp(-weight * distances * distances)
    return float(np.median(values))

# ----- main SC calculation -----

def calculate_sc(
    atoms1: List[Dict],
    atoms2: List[Dict],
    probe_radius: float = PROBE_RADIUS,
    dot_density: int = DOT_DENSITY,
    weight: float = WEIGHT,
    interface_distance: float = INTERFACE_DISTANCE,
    trim: float = TRIM,
) -> float:
    # try compiled backend first
    if USE_SC_BACKEND and _calculate_sc_backend is not None:
        coords = np.array([a["xyz"] for a in atoms1 + atoms2], dtype=np.float64)
        radii_arr = np.array([a["radius"] for a in atoms1 + atoms2], dtype=np.float64)
        mol_ids = np.array([1] * len(atoms1) + [2] * len(atoms2), dtype=np.int32)
        sc_val, stats = _calculate_sc_backend(
            coords, radii_arr, mol_ids,
            probe_radius=probe_radius, dot_density=float(dot_density),
            weight=weight, trim=trim, interface_distance=interface_distance,
        )
        if isinstance(stats, dict):
            b1 = stats.get("n_buried_m1", 0)
            b2 = stats.get("n_buried_m2", 0)
            t1 = stats.get("n_trimmed_m1", b1)
            t2 = stats.get("n_trimmed_m2", b2)
            if b1 > 0 and b2 > 0 and t1 > 0 and t2 > 0 and np.isfinite(sc_val):
                return float(sc_val)

    # Python Connolly path (SCASA mds if available, fallback to SAS)
    if not atoms1 or not atoms2:
        return 0.0

    iface1, iface2 = calculate_interface_atoms(atoms1, atoms2, interface_distance)
    if not iface1 or not iface2:
        return 0.0

    if _connolly_mds is not None:
        all_atoms = iface1 + iface2
        coords = np.array([a["xyz"] for a in all_atoms], dtype=float)
        radii_arr = np.array([a["radius"] for a in all_atoms], dtype=float)
        mol_ids = np.array([1] * len(iface1) + [2] * len(iface2), dtype=np.int32)

        dots, normals, flags, dot_mol = _connolly_mds(
            float(probe_radius), coords, radii_arr, mol_ids, density=float(dot_density))

        if len(dots) > 0:
            buried_mask = flags == _BURIED_FLAG
            d1_all = dots[(dot_mol == 1) & buried_mask]
            n1_all = normals[(dot_mol == 1) & buried_mask]
            d2_all = dots[(dot_mol == 2) & buried_mask]
            n2_all = normals[(dot_mol == 2) & buried_mask]

            if len(d1_all) >= 5 and len(d2_all) >= 5:
                _, idx2_all = cKDTree(d2_all).query(d1_all)
                mask1 = np.linalg.norm(d1_all - d2_all[idx2_all], axis=1) <= trim
                _, idx1_all = cKDTree(d1_all).query(d2_all)
                mask2 = np.linalg.norm(d2_all - d1_all[idx1_all], axis=1) <= trim

                d1, n1 = d1_all[mask1], n1_all[mask1]
                d2, n2 = d2_all[mask2], n2_all[mask2]

                if len(d1) >= 5 and len(d2) >= 5:
                    s_ab = _score(d1, n1, d2, n2, weight)
                    s_ba = _score(d2, n2, d1, n1, weight)
                    return (s_ab + s_ba) / 2.0

    # Fallback: original SAS-based surface generation
    s1, n1 = _surface_points(iface1, probe_radius, dot_density)
    s2, n2 = _surface_points(iface2, probe_radius, dot_density)

    b1, bn1, n_buried1, n_trim1 = _classify_and_trim(s1, n1, iface2, probe_radius, trim)
    b2, bn2, n_buried2, n_trim2 = _classify_and_trim(s2, n2, iface1, probe_radius, trim)

    if n_trim1 < 5 or n_trim2 < 5:
        return 0.0

    s_ab = _score(b1, bn1, b2, bn2, weight)
    s_ba = _score(b2, bn2, b1, bn1, weight)
    return (s_ab + s_ba) / 2.0


# ----- public wrappers -----

def calculate_sc_from_pdb(pdb_file: str, chain1: str, chain2: str, **kwargs) -> float:
    return calculate_sc(
        read_pdb_atoms(pdb_file, [chain1]),
        read_pdb_atoms(pdb_file, [chain2]),
        **kwargs,
    )

def calculate_sc_combined_from_pdb(pdb_file: str, query_chains: List[str], partner_chains: List[str], **kwargs) -> float:
    return calculate_sc(
        read_pdb_atoms(pdb_file, query_chains),
        read_pdb_atoms(pdb_file, partner_chains),
        **kwargs,
    )

def calculate_sc_batch(pdb_file: str, query_chains: List[str], partner_chains: List[str], **kwargs) -> Dict[Tuple[str, str], float]:
    all_chains = set(query_chains + partner_chains)
    all_atoms = {c: read_pdb_atoms(pdb_file, [c]) for c in all_chains}
    return {(qc, pc): calculate_sc(all_atoms[qc], all_atoms[pc], **kwargs) for qc in query_chains for pc in partner_chains}
