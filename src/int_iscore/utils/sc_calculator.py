"""Shape Complementarity (SC) calculator rebuilt around a voxel SES model."""

from __future__ import annotations

from fnmatch import fnmatchcase
from functools import lru_cache
import os
from typing import Dict, List, Tuple

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

try:
    from .sc_backend import calculate_sc_backend as _calculate_sc_backend
except Exception:
    _calculate_sc_backend = None

USE_SC_BACKEND = os.getenv("INT_ISCORE_USE_SC_BACKEND", "0") == "1"


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
GRID_SPACING = 1.0 / np.sqrt(DOT_DENSITY)
PADDING = 3.0
MAX_GRID_POINTS = 24_000_000
DIST_CHUNK = 512


def _parse_sc_radii() -> List[Tuple[str, str, float]]:
    rules = []
    for line in SC_RADII_LIB.splitlines():
        parts = line.split()
        if len(parts) == 3:
            rules.append((parts[0], parts[1], float(parts[2])))
    return rules


RADII_RULES = _parse_sc_radii()


def get_atom_radius(atom_name: str, res_name: str) -> float:
    atom_name = atom_name.strip()
    res_name = res_name.strip()
    for res_pat, atom_pat, radius in RADII_RULES:
        if fnmatchcase(res_name, res_pat) and fnmatchcase(atom_name, atom_pat):
            return radius
    return DEFAULT_RADIUS


def read_pdb_atoms(pdb_file: str, chains: List[str]) -> List[Dict]:
    atoms = []
    with open(pdb_file, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = line[21:22].strip()
            if chain not in chains:
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            atoms.append(
                {
                    "name": atom_name,
                    "res": line[17:20].strip(),
                    "res_id": int(line[22:26].strip()),
                    "chain": chain,
                    "xyz": np.array(
                        [
                            float(line[30:38].strip()),
                            float(line[38:46].strip()),
                            float(line[46:54].strip()),
                        ],
                        dtype=float,
                    ),
                    "radius": get_atom_radius(atom_name, line[17:20].strip()),
                }
            )
    return atoms


@lru_cache(maxsize=32)
def _sphere_structure(radius_voxels: int) -> np.ndarray:
    coords = np.arange(-radius_voxels, radius_voxels + 1)
    x, y, z = np.meshgrid(coords, coords, coords, indexing="ij")
    return (x * x + y * y + z * z) <= (radius_voxels * radius_voxels)


def calculate_interface_atoms(atoms1: List[Dict], atoms2: List[Dict], interface_distance: float = INTERFACE_DISTANCE) -> Tuple[List[Dict], List[Dict]]:
    if not atoms1 or not atoms2:
        return [], []
    coords1 = np.array([a["xyz"] for a in atoms1], dtype=float)
    coords2 = np.array([a["xyz"] for a in atoms2], dtype=float)
    radii1 = np.array([a["radius"] for a in atoms1], dtype=float)
    radii2 = np.array([a["radius"] for a in atoms2], dtype=float)
    distances = np.linalg.norm(coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :], axis=2)
    clearance = distances - radii1[:, np.newaxis] - radii2[np.newaxis, :]
    keep1 = np.min(clearance, axis=1) < interface_distance
    keep2 = np.min(clearance, axis=0) < interface_distance
    return [atom for atom, keep in zip(atoms1, keep1) if keep], [atom for atom, keep in zip(atoms2, keep2) if keep]


def _make_vdw_volume(atoms: List[Dict], spacing: float, padding: float) -> Tuple[np.ndarray, np.ndarray]:
    coords = np.array([atom["xyz"] for atom in atoms], dtype=float)
    radii = np.array([atom["radius"] for atom in atoms], dtype=float)
    min_corner = np.min(coords - radii[:, np.newaxis], axis=0) - padding
    max_corner = np.max(coords + radii[:, np.newaxis], axis=0) + padding
    shape = np.ceil((max_corner - min_corner) / spacing).astype(int) + 3
    volume = np.zeros(tuple(shape.tolist()), dtype=bool)
    for atom in atoms:
        r = atom["radius"]
        r_vox = int(np.ceil(r / spacing))
        center = ((atom["xyz"] - min_corner) / spacing).astype(int)
        x0, y0, z0 = np.maximum(center - r_vox - 1, 0)
        x1, y1, z1 = np.minimum(center + r_vox + 2, shape)
        xs = np.arange(x0, x1)
        ys = np.arange(y0, y1)
        zs = np.arange(z0, z1)
        xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
        points = min_corner + np.stack((xx, yy, zz), axis=-1) * spacing
        mask = np.sum((points - atom["xyz"]) ** 2, axis=-1) <= r * r
        volume[x0:x1, y0:y1, z0:z1] |= mask
    return volume, min_corner


def _choose_grid_spacing(atoms: List[Dict], probe_radius: float, target_spacing: float) -> float:
    coords = np.array([atom["xyz"] for atom in atoms], dtype=float)
    radii = np.array([atom["radius"] for atom in atoms], dtype=float)
    min_corner = np.min(coords - radii[:, np.newaxis], axis=0) - (PADDING + probe_radius)
    max_corner = np.max(coords + radii[:, np.newaxis], axis=0) + (PADDING + probe_radius)
    extent = np.maximum(max_corner - min_corner, 1e-6)
    spacing = target_spacing
    est_points = np.prod(np.ceil(extent / spacing).astype(int) + 3)
    while est_points > MAX_GRID_POINTS:
        spacing *= 1.15
        est_points = np.prod(np.ceil(extent / spacing).astype(int) + 3)
    return float(spacing)


def _build_ses_surface(atoms: List[Dict], probe_radius: float, spacing: float) -> Tuple[np.ndarray, np.ndarray]:
    if not atoms:
        return np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=float)
    vdw_volume, origin = _make_vdw_volume(atoms, spacing, PADDING + probe_radius)
    probe_vox = max(1, int(np.ceil(probe_radius / spacing)))
    closed = ndi.binary_closing(vdw_volume, structure=_sphere_structure(probe_vox))
    boundary = closed & ~ndi.binary_erosion(closed, structure=np.ones((3, 3, 3), dtype=bool))
    if not np.any(boundary):
        return np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=float)

    inside = ndi.distance_transform_edt(closed, sampling=spacing)
    outside = ndi.distance_transform_edt(~closed, sampling=spacing)
    signed = inside - outside
    grad = np.stack(np.gradient(signed, spacing), axis=-1)

    idx = np.argwhere(boundary)
    points = origin + idx * spacing
    normals = grad[idx[:, 0], idx[:, 1], idx[:, 2]]
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / np.clip(norms, 1e-8, None)
    return points, normals


def _build_surface_with_auto_spacing(atoms: List[Dict], probe_radius: float) -> Tuple[np.ndarray, np.ndarray, float]:
    spacing = _choose_grid_spacing(atoms, probe_radius, GRID_SPACING)
    points, normals = _build_ses_surface(atoms, probe_radius, spacing)
    return points, normals, spacing


def _classify_surface_points(points: np.ndarray, normals: np.ndarray, complex_surface_points: np.ndarray, threshold: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if len(points) == 0:
        empty = np.empty((0, 3), dtype=float)
        return empty, empty, empty, empty
    if len(complex_surface_points) == 0:
        empty = np.empty((0, 3), dtype=float)
        return empty, empty, points, normals
    mins = _min_distances_chunked(points, complex_surface_points)
    buried = mins > threshold
    return points[buried], normals[buried], points[~buried], normals[~buried]


def _min_distances_chunked(points_a: np.ndarray, points_b: np.ndarray) -> np.ndarray:
    if len(points_a) == 0:
        return np.empty((0,), dtype=float)
    if len(points_b) == 0:
        return np.full(len(points_a), np.inf, dtype=float)
    tree = cKDTree(points_b)
    mins, _ = tree.query(points_a, k=1, workers=-1)
    return mins.astype(float, copy=False)


def _trim_peripheral_band(buried_points: np.ndarray, buried_normals: np.ndarray, accessible_points: np.ndarray, trim: float) -> Tuple[np.ndarray, np.ndarray]:
    if len(buried_points) == 0 or len(accessible_points) == 0:
        return buried_points, buried_normals
    min_distances = _min_distances_chunked(buried_points, accessible_points)
    keep = min_distances >= trim
    if np.count_nonzero(keep) >= 64:
        return buried_points[keep], buried_normals[keep]

    # Thin interfaces can be over-trimmed by a fixed peripheral band. Relax the
    # threshold conservatively until we keep a minimally stable point set.
    relaxed_thresholds = [trim * 0.8, trim * 0.6, trim * 0.4, trim * 0.2]
    for threshold in relaxed_thresholds:
        keep = min_distances >= threshold
        if np.count_nonzero(keep) >= 64:
            return buried_points[keep], buried_normals[keep]

    # Final fallback: keep the most interior buried points rather than returning
    # an empty interface and collapsing SC to zero.
    keep_count = min(len(buried_points), 256)
    best = np.argsort(min_distances)[-keep_count:]
    return buried_points[best], buried_normals[best]


def _classify_and_trim(points: np.ndarray, normals: np.ndarray, complex_surface_points: np.ndarray, threshold: float, trim: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    buried_points, buried_normals, accessible_points, accessible_normals = _classify_surface_points(points, normals, complex_surface_points, threshold)
    trimmed_points, trimmed_normals = _trim_peripheral_band(buried_points, buried_normals, accessible_points, trim)
    return buried_points, buried_normals, accessible_points, accessible_normals, trimmed_points, trimmed_normals


def _surface_complementarity_values(interface_points: np.ndarray, interface_normals: np.ndarray, full_points: np.ndarray, full_normals: np.ndarray, weight: float) -> List[float]:
    if len(interface_points) == 0 or len(full_points) == 0:
        return []
    tree = cKDTree(full_points)
    min_distances, nearest = tree.query(interface_points, k=1, workers=-1)
    dot = -np.sum(interface_normals * full_normals[nearest], axis=1)
    return (dot * np.exp(-weight * min_distances * min_distances)).tolist()


def calculate_sc(
    atoms1: List[Dict],
    atoms2: List[Dict],
    probe_radius: float = PROBE_RADIUS,
    dot_density: int = DOT_DENSITY,
    weight: float = WEIGHT,
    interface_distance: float = INTERFACE_DISTANCE,
    trim: float = TRIM,
) -> float:
    if USE_SC_BACKEND and _calculate_sc_backend is not None:
        coords = np.array([a["xyz"] for a in atoms1 + atoms2], dtype=np.float64)
        radii = np.array([a["radius"] for a in atoms1 + atoms2], dtype=np.float64)
        molecule_ids = np.array([1] * len(atoms1) + [2] * len(atoms2), dtype=np.int32)
        sc_value, _stats = _calculate_sc_backend(
            coords,
            radii,
            molecule_ids,
            probe_radius=probe_radius,
            dot_density=float(dot_density),
            weight=weight,
            trim=trim,
            interface_distance=interface_distance,
        )
        if isinstance(_stats, dict):
            buried1 = int(_stats.get("n_buried_m1", 0))
            buried2 = int(_stats.get("n_buried_m2", 0))
            trimmed1 = int(_stats.get("n_trimmed_m1", buried1))
            trimmed2 = int(_stats.get("n_trimmed_m2", buried2))
            if buried1 > 0 and buried2 > 0 and trimmed1 > 0 and trimmed2 > 0 and np.isfinite(sc_value):
                return float(sc_value)

        # Backend is experimental; fall back to the current Python reference path
        # when it fails to produce a plausible interfacial surface.

    del dot_density
    interface_atoms1, interface_atoms2 = calculate_interface_atoms(atoms1, atoms2, interface_distance)
    if not interface_atoms1 or not interface_atoms2:
        return 0.0
    combined_atoms = interface_atoms1 + interface_atoms2
    spacing = _choose_grid_spacing(combined_atoms, probe_radius, GRID_SPACING)
    surface1, normals1 = _build_ses_surface(interface_atoms1, probe_radius, spacing)
    surface2, normals2 = _build_ses_surface(interface_atoms2, probe_radius, spacing)
    complex_surface, _ = _build_ses_surface(combined_atoms, probe_radius, spacing)

    threshold = spacing * 1.25
    buried1, buried_normals1, accessible1, _, trimmed1, trimmed_normals1 = _classify_and_trim(surface1, normals1, complex_surface, threshold, trim)
    buried2, buried_normals2, accessible2, _, trimmed2, trimmed_normals2 = _classify_and_trim(surface2, normals2, complex_surface, threshold, trim)

    if len(trimmed1) >= 64:
        buried1, buried_normals1 = trimmed1, trimmed_normals1
    if len(trimmed2) >= 64:
        buried2, buried_normals2 = trimmed2, trimmed_normals2

    if len(buried1) < 64 or len(buried2) < 64:
        for factor in (1.0, 1.5, 2.0, 2.5):
            alt_b1, alt_bn1, alt_a1, _, alt_t1, alt_tn1 = _classify_and_trim(surface1, normals1, complex_surface, spacing * factor, trim)
            alt_b2, alt_bn2, alt_a2, _, alt_t2, alt_tn2 = _classify_and_trim(surface2, normals2, complex_surface, spacing * factor, trim)
            cand1 = alt_t1 if len(alt_t1) >= 64 else alt_b1
            candn1 = alt_tn1 if len(alt_t1) >= 64 else alt_bn1
            cand2 = alt_t2 if len(alt_t2) >= 64 else alt_b2
            candn2 = alt_tn2 if len(alt_t2) >= 64 else alt_bn2
            if len(cand1) > len(buried1):
                buried1, buried_normals1 = cand1, candn1
            if len(cand2) > len(buried2):
                buried2, buried_normals2 = cand2, candn2
    if len(buried1) < 5 or len(buried2) < 5:
        return 0.0
    s_ab = _surface_complementarity_values(buried1, buried_normals1, surface2, normals2, weight)
    s_ba = _surface_complementarity_values(buried2, buried_normals2, surface1, normals1, weight)
    if len(s_ab) < 5 or len(s_ba) < 5:
        return 0.0
    return (float(np.median(s_ab)) + float(np.median(s_ba))) / 2.0


def calculate_sc_from_pdb(pdb_file: str, chain1: str, chain2: str, **kwargs) -> float:
    return calculate_sc(read_pdb_atoms(pdb_file, [chain1]), read_pdb_atoms(pdb_file, [chain2]), **kwargs)


def calculate_sc_combined_from_pdb(pdb_file: str, query_chains: List[str], partner_chains: List[str], **kwargs) -> float:
    return calculate_sc(read_pdb_atoms(pdb_file, query_chains), read_pdb_atoms(pdb_file, partner_chains), **kwargs)


def calculate_sc_batch(pdb_file: str, query_chains: List[str], partner_chains: List[str], **kwargs) -> Dict[Tuple[str, str], float]:
    all_chains = set(query_chains + partner_chains)
    all_atoms = {chain: read_pdb_atoms(pdb_file, [chain]) for chain in all_chains}
    return {(qc, pc): calculate_sc(all_atoms[qc], all_atoms[pc], **kwargs) for qc in query_chains for pc in partner_chains}


def calculate_sc_summary(pdb_file: str, query_chains: List[str], partner_chains: List[str], **kwargs) -> Dict:
    sc_scores = calculate_sc_batch(pdb_file, query_chains, partner_chains, **kwargs)
    values = list(sc_scores.values()) if sc_scores else [0.0]
    return {
        "sc_scores": sc_scores,
        "n_pairs": len(sc_scores),
        "sc_mean": float(np.mean(values)) if sc_scores else 0.0,
        "sc_min": float(np.min(values)) if sc_scores else 0.0,
        "sc_max": float(np.max(values)) if sc_scores else 0.0,
    }
