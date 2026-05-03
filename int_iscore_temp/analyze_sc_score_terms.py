from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, "/media/cuixi/data01/lluoto/int-iScore/src")
sys.path.insert(0, "/media/cuixi/data01/lluoto/int-iScore")

from int_iscore.utils import sc_calculator as sc


def inspect(pdb_path: str):
    at1 = sc.read_pdb_atoms(pdb_path, ["A"])
    at2 = sc.read_pdb_atoms(pdb_path, ["B"])
    i1, i2 = sc.calculate_interface_atoms(at1, at2, sc.INTERFACE_DISTANCE)
    combined = i1 + i2
    spacing = sc._choose_grid_spacing(combined, sc.PROBE_RADIUS, sc.GRID_SPACING)
    s1, n1 = sc._build_ses_surface(i1, sc.PROBE_RADIUS, spacing)
    s2, n2 = sc._build_ses_surface(i2, sc.PROBE_RADIUS, spacing)
    scx, _ = sc._build_ses_surface(combined, sc.PROBE_RADIUS, spacing)
    thr = spacing * 1.25
    b1, bn1, a1, _ = sc._classify_surface_points(s1, n1, scx, thr)
    b2, bn2, a2, _ = sc._classify_surface_points(s2, n2, scx, thr)
    t1, tn1 = sc._trim_peripheral_band(b1, bn1, a1, sc.TRIM)
    t2, tn2 = sc._trim_peripheral_band(b2, bn2, a2, sc.TRIM)
    if len(t1) < 64:
        t1, tn1 = b1, bn1
    if len(t2) < 64:
        t2, tn2 = b2, bn2

    from scipy.spatial import cKDTree
    tree2 = cKDTree(s2)
    d12, idx12 = tree2.query(t1, k=1, workers=-1)
    dot12 = -np.sum(tn1 * n2[idx12], axis=1)
    s12 = dot12 * np.exp(-sc.WEIGHT * d12 * d12)

    tree1 = cKDTree(s1)
    d21, idx21 = tree1.query(t2, k=1, workers=-1)
    dot21 = -np.sum(tn2 * n1[idx21], axis=1)
    s21 = dot21 * np.exp(-sc.WEIGHT * d21 * d21)

    print(pdb_path)
    print('trim', len(t1), len(t2))
    print('d12 mean/med', float(np.mean(d12)), float(np.median(d12)))
    print('d21 mean/med', float(np.mean(d21)), float(np.median(d21)))
    print('dot12 mean/med', float(np.mean(dot12)), float(np.median(dot12)))
    print('dot21 mean/med', float(np.mean(dot21)), float(np.median(dot21)))
    print('s12 mean/med', float(np.mean(s12)), float(np.median(s12)))
    print('s21 mean/med', float(np.mean(s21)), float(np.median(s21)))
    print('sc', float((np.median(s12) + np.median(s21)) / 2.0))
    print()


for path in [
    '/media/cuixi/data01/lluoto/6A6I/seed-6893703_sample-1/model.pdb',
    '/media/cuixi/data01/lluoto/6A6I/seed-958038_sample-1/model.pdb',
    '/media/cuixi/data01/lluoto/6A6I/seed-614460_sample-1/model.pdb',
    '/media/cuixi/data01/lluoto/5-HT2C-2CD28/pred.rank_0.pdb',
]:
    inspect(path)
