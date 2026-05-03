import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore/src")
sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore")

from int_iscore.utils import sc_calculator as sc


HELPER = "/mnt/sdb2/lluoto/int-iScore/int_iscore_temp/run_ccp4_sc.sh"


def ccp4_stats(pdb: str, chain1: str, *chains2: str):
    proc = subprocess.run(
        ["bash", HELPER, pdb, chain1, *chains2],
        check=True,
        capture_output=True,
        text=True,
    )
    text = proc.stdout + "\n" + proc.stderr
    stats = {}
    patterns = {
        "potential_interface_atoms": r"Potential interface atoms\s+(\d+)",
        "total_surface_points": r"Total\s+(\d+)\s*$",
        "dots_mol1": r"Number of dots read for molecule 1 =\s+(\d+)",
        "dots_mol2": r"Number of dots read for molecule 2 =\s+(\d+)",
        "buried_mol1": r"Number of dots buried for molecule 1 =\s+(\d+)",
        "buried_mol2": r"Number of dots buried for molecule 2 =\s+(\d+)",
        "accessible_mol1": r"Number of dots accessible for molecule 1 =\s+(\d+)",
        "accessible_mol2": r"Number of dots accessible for molecule 2 =\s+(\d+)",
        "trimmed_mol1": r"Number of points left for molecule 1 after trim is\s+(\d+)",
        "trimmed_mol2": r"Number of points left for molecule 2 after trim is\s+(\d+)",
        "sc": r"Shape complementarity statistic Sc =\s+([0-9.\-]+)",
    }
    for key, pattern in patterns.items():
        m = re.search(pattern, text, re.MULTILINE)
        if m:
            stats[key] = float(m.group(1)) if key == "sc" else int(m.group(1))
        else:
            stats[key] = None
    return stats


def py_stats(pdb: str, chains1, chains2):
    at1 = sc.read_pdb_atoms(pdb, list(chains1))
    at2 = sc.read_pdb_atoms(pdb, list(chains2))
    i1, i2 = sc.calculate_interface_atoms(at1, at2, sc.INTERFACE_DISTANCE)
    combined_atoms = i1 + i2
    spacing = sc._choose_grid_spacing(combined_atoms, sc.PROBE_RADIUS, sc.GRID_SPACING)
    s1, n1 = sc._build_ses_surface(i1, sc.PROBE_RADIUS, spacing)
    s2, n2 = sc._build_ses_surface(i2, sc.PROBE_RADIUS, spacing)
    b1, bn1, a1, _ = sc._classify_surface_points(s1, n1, i2, sc.PROBE_RADIUS)
    b2, bn2, a2, _ = sc._classify_surface_points(s2, n2, i1, sc.PROBE_RADIUS)
    t1, tn1 = sc._trim_peripheral_band(b1, bn1, a1, sc.TRIM)
    t2, tn2 = sc._trim_peripheral_band(b2, bn2, a2, sc.TRIM)
    if len(t1) < 64 or len(t2) < 64:
        complex_surface, _ = sc._build_ses_surface(combined_atoms, sc.PROBE_RADIUS, spacing)
        alt1, alt_normals1 = sc._classify_by_complex_surface(s1, n1, complex_surface, spacing * 2.0, sc.TRIM)
        alt2, alt_normals2 = sc._classify_by_complex_surface(s2, n2, complex_surface, spacing * 2.0, sc.TRIM)
        if len(alt1) >= len(t1):
            t1, tn1 = alt1, alt_normals1
        if len(alt2) >= len(t2):
            t2, tn2 = alt2, alt_normals2
    return {
        "interface_atoms_mol1": len(i1),
        "interface_atoms_mol2": len(i2),
        "spacing": spacing,
        "surface_mol1": len(s1),
        "surface_mol2": len(s2),
        "buried_mol1": len(b1),
        "buried_mol2": len(b2),
        "accessible_mol1": len(a1),
        "accessible_mol2": len(a2),
        "trimmed_mol1": len(t1),
        "trimmed_mol2": len(t2),
        "sc": sc.calculate_sc(at1, at2),
    }


cases = [
    ("5HT2C A-B", "/mnt/sdb2/lluoto/5-HT2C-2CD28/pred.rank_0.pdb", ["A"], ["B"]),
    ("6A6I outlier A-B", "/mnt/sdb2/lluoto/6A6I/seed-9045222_sample-0/model.pdb", ["A"], ["B"]),
]

for title, pdb, q, p in cases:
    print(f"## {title}")
    print("CCP4")
    print(ccp4_stats(pdb, q[0], *p))
    print("Python")
    print(py_stats(pdb, q, p))
    print()
