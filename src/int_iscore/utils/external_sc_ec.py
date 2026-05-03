"""SC/EC helpers without CCP4 or Chimera dependencies."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np

from .sc_calculator import calculate_sc_batch, calculate_sc_from_pdb

PROBE_RADIUS = 1.4
SHELL_POINTS = 96
SHELL_TOLERANCE = 1e-3


def _resolve_binary(explicit: Optional[str], candidates: Iterable[str]) -> Optional[str]:
    if explicit:
        explicit_path = Path(explicit)
        if explicit_path.is_file():
            return str(explicit_path)
        found = shutil.which(explicit)
        if found:
            return found
    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return found
    return None


def _resolve_apbs_binary(explicit: Optional[str]) -> Optional[str]:
    resolved = _resolve_binary(explicit, ("apbs",))
    if resolved:
        return resolved
    fallback = Path("/home/cuixi/miniforge3/bin/apbs")
    if fallback.is_file():
        return str(fallback)
    return None


def _resolve_env_binary(explicit: Optional[str], name: str) -> Optional[str]:
    resolved = _resolve_binary(explicit, (name,))
    if resolved:
        return resolved
    fallback = Path(f"/home/cuixi/miniforge3/envs/int_iScore/bin/{name}")
    if fallback.is_file():
        return str(fallback)
    return None


def calculate_sc_external(
    pdb_file: str,
    query_chains,
    partner_chains,
    ccp4_path: Optional[str] = None,
) -> Optional[float]:
    """Compatibility wrapper that now uses the in-repo SC implementation."""
    del ccp4_path
    if not query_chains or not partner_chains:
        return None
    if len(query_chains) == 1 and len(partner_chains) == 1:
        return calculate_sc_from_pdb(pdb_file, query_chains[0], partner_chains[0])

    sc_scores = calculate_sc_batch(pdb_file, query_chains, partner_chains)
    if not sc_scores:
        return None
    return float(np.mean(list(sc_scores.values())))


def _read_pqr_atoms(pqr_file: str):
    atoms = []
    with open(pqr_file, "r") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            atoms.append(
                {
                    "atom_name": atom_name,
                    "res_id": int(line[22:26].strip()),
                    "chain": line[21:22].strip(),
                    "xyz": np.array(
                        [
                            float(line[30:38].strip()),
                            float(line[38:46].strip()),
                            float(line[46:54].strip()),
                        ],
                        dtype=float,
                    ),
                    "radius": float(line[62:70].strip()),
                }
            )
    return atoms


def _parse_dx_grid(dx_file: str):
    counts = None
    origin = None
    deltas = []
    values = []
    reading_values = False

    with open(dx_file, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("object 1 class gridpositions counts"):
                counts = tuple(int(x) for x in line.split()[-3:])
                continue
            if line.startswith("origin"):
                origin = np.array([float(x) for x in line.split()[1:4]], dtype=float)
                continue
            if line.startswith("delta"):
                deltas.append(np.array([float(x) for x in line.split()[1:4]], dtype=float))
                continue
            if "data follows" in line:
                reading_values = True
                continue
            if reading_values:
                if line.startswith("attribute") or line.startswith("object"):
                    reading_values = False
                    continue
                values.extend(float(x) for x in line.split())

    if counts is None or origin is None or len(deltas) != 3:
        raise ValueError(f"Incomplete DX grid header in {dx_file}")

    axes = []
    spacing = []
    for delta in deltas:
        nz = np.flatnonzero(np.abs(delta) > 1e-12)
        if len(nz) != 1:
            raise ValueError("Only axis-aligned APBS DX grids are supported")
        axes.append(int(nz[0]))
        spacing.append(float(delta[nz[0]]))

    if axes != [0, 1, 2]:
        raise ValueError("Unexpected DX axis order")

    data = np.asarray(values, dtype=float)
    expected = int(np.prod(counts))
    if data.size != expected:
        raise ValueError(f"DX grid has {data.size} values, expected {expected}")

    return {
        "origin": origin,
        "counts": np.array(counts, dtype=int),
        "spacing": np.array(spacing, dtype=float),
        "data": data.reshape(counts, order="C"),
    }


def _fibonacci_sphere_points(n_points: int) -> np.ndarray:
    if n_points <= 1:
        return np.array([[0.0, 0.0, 1.0]], dtype=float)

    indices = np.arange(n_points, dtype=float)
    phi = np.pi * (3.0 - np.sqrt(5.0))
    y = 1.0 - (2.0 * indices) / (n_points - 1)
    radius = np.sqrt(np.clip(1.0 - y * y, 0.0, None))
    theta = phi * indices
    x = np.cos(theta) * radius
    z = np.sin(theta) * radius
    return np.column_stack((x, y, z))


UNIT_SPHERE_POINTS = _fibonacci_sphere_points(SHELL_POINTS)


def _interpolate_dx_values(points: np.ndarray, grid: Dict) -> np.ndarray:
    origin = grid["origin"]
    spacing = grid["spacing"]
    counts = grid["counts"]
    data = grid["data"]

    fractional = (points - origin) / spacing
    inside = np.all((fractional >= 0.0) & (fractional < (counts - 1)), axis=1)
    values = np.full(len(points), np.nan, dtype=float)
    if not np.any(inside):
        return values

    frac = fractional[inside]
    base = np.floor(frac).astype(int)
    delta = frac - base
    x0, y0, z0 = base.T
    x1 = x0 + 1
    y1 = y0 + 1
    z1 = z0 + 1
    dx, dy, dz = delta.T

    c000 = data[x0, y0, z0]
    c100 = data[x1, y0, z0]
    c010 = data[x0, y1, z0]
    c110 = data[x1, y1, z0]
    c001 = data[x0, y0, z1]
    c101 = data[x1, y0, z1]
    c011 = data[x0, y1, z1]
    c111 = data[x1, y1, z1]

    c00 = c000 * (1.0 - dx) + c100 * dx
    c10 = c010 * (1.0 - dx) + c110 * dx
    c01 = c001 * (1.0 - dx) + c101 * dx
    c11 = c011 * (1.0 - dx) + c111 * dx
    c0 = c00 * (1.0 - dy) + c10 * dy
    c1 = c01 * (1.0 - dy) + c11 * dy
    values[inside] = c0 * (1.0 - dz) + c1 * dz
    return values


def _accessible_shell_potentials(atom: Dict, all_coords: np.ndarray, all_radii: np.ndarray, grid: Dict) -> np.ndarray:
    shell_radius = atom["radius"] + PROBE_RADIUS
    shell_points = atom["xyz"] + UNIT_SPHERE_POINTS * shell_radius

    distances = np.linalg.norm(shell_points[:, np.newaxis, :] - all_coords[np.newaxis, :, :], axis=2)
    accessible = np.all(
        distances >= (all_radii[np.newaxis, :] + PROBE_RADIUS - SHELL_TOLERANCE),
        axis=1,
    )
    accessible_points = shell_points[accessible]
    if len(accessible_points) == 0:
        return np.array([], dtype=float)

    sampled = _interpolate_dx_values(accessible_points, grid)
    sampled = sampled[np.isfinite(sampled)]
    return sampled


def _calculate_ec_score(residue_potentials_dict: Dict[str, list[float]]) -> float:
    all_numerator = 0.0
    all_denominator = 0.0
    for resi_own in residue_potentials_dict.values():
        if len(resi_own) <= 1:
            continue
        phi_array = np.asarray(resi_own, dtype=float)
        phi_mean = float(np.mean(phi_array))
        for index, phi_current in enumerate(phi_array):
            phi_rest = np.delete(phi_array, index)
            phi_prime_mean = float(np.mean(phi_rest))
            numerator = (phi_current - phi_mean) * (np.sum(phi_rest) - phi_prime_mean)
            denominator = np.sqrt(
                np.sum((phi_current - phi_mean) ** 2)
                * np.sum((phi_rest - phi_prime_mean) ** 2)
            )
            all_numerator += float(numerator)
            all_denominator += float(denominator)
    return float(-all_numerator / all_denominator) if all_denominator != 0 else 0.0


def _map_atoms_python(pqr_file: str, interface_residues, dx_file: str) -> Optional[float]:
    interface_residue_set = set(interface_residues)
    residue_potentials_dict = {residue: [] for residue in interface_residue_set}
    pqr_atoms = _read_pqr_atoms(pqr_file)
    if not pqr_atoms:
        return None

    grid = _parse_dx_grid(dx_file)
    all_coords = np.array([atom["xyz"] for atom in pqr_atoms], dtype=float)
    all_radii = np.array([atom["radius"] for atom in pqr_atoms], dtype=float)

    for atom in pqr_atoms:
        residue_key = f'{atom["res_id"]}_{atom["chain"]}'
        if residue_key not in residue_potentials_dict:
            continue
        sampled = _accessible_shell_potentials(atom, all_coords, all_radii, grid)
        if sampled.size:
            residue_potentials_dict[residue_key].extend(sampled.tolist())

    if not any(values for values in residue_potentials_dict.values()):
        return None
    return _calculate_ec_score(residue_potentials_dict)


def calculate_ec_external(
    pdb_file: str,
    interface_residues,
    chimera_path: Optional[str] = None,
    vmd_path: Optional[str] = None,
    pdb2pqr_path: Optional[str] = None,
    apbs_path: Optional[str] = None,
) -> Optional[float]:
    """Calculate EC from APBS output using direct Python grid sampling."""
    del chimera_path, vmd_path
    if not interface_residues:
        return None

    pdbfixer_bin = _resolve_env_binary(None, "pdbfixer")
    pdb2pqr_bin = _resolve_binary(pdb2pqr_path, ("pdb2pqr", "pdb2pqr30")) or _resolve_env_binary(None, "pdb2pqr")
    apbs_bin = _resolve_apbs_binary(apbs_path)
    if not pdbfixer_bin or not pdb2pqr_bin or not apbs_bin:
        return None

    with tempfile.TemporaryDirectory(prefix="intiscore_ec_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        pdb_src = Path(pdb_file)
        pdb_copy = tmpdir_path / pdb_src.name
        shutil.copy2(pdb_src, pdb_copy)
        refined = tmpdir_path / f"{pdb_src.stem}_refine.pdb"
        pqr_file = tmpdir_path / f"{pdb_src.stem}.pqr"
        in_file = tmpdir_path / f"{pdb_src.stem}.in"
        dx_file = Path(str(pqr_file) + ".dx")

        subprocess.run(
            [pdbfixer_bin, str(pdb_copy), "--output", str(refined)],
            check=True,
            capture_output=True,
            text=True,
            cwd=tmpdir_path,
        )
        subprocess.run(
            [
                pdb2pqr_bin,
                "--noopt",
                "--ff",
                "AMBER",
                "--keep-chain",
                "--apbs-input",
                str(in_file),
                str(refined),
                str(pqr_file),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=tmpdir_path,
        )
        subprocess.run([apbs_bin, str(in_file)], check=True, capture_output=True, text=True, cwd=tmpdir_path)

        return _map_atoms_python(str(pqr_file), interface_residues, str(dx_file))
