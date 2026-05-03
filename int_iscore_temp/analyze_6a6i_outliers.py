import math
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, "/media/cuixi/data01/lluoto/int-iScore/src")
sys.path.insert(0, "/media/cuixi/data01/lluoto/int-iScore")

from int_iscore.utils import sc_calculator as sc


BASE = Path("/media/cuixi/data01/lluoto/6A6I")
CSV = Path("/media/cuixi/data01/lluoto/int-iScore/int_iscore_temp/6A6I_sc_comparison.csv")


def filename_to_model_path(name: str) -> Path:
    suffix = name[len("6A6I_"):]
    stem = suffix[:-len("model")]
    return BASE / stem / "model.pdb"


def surface_stats(pdb: Path):
    at1 = sc.read_pdb_atoms(str(pdb), ["A"])
    at2 = sc.read_pdb_atoms(str(pdb), ["B"])
    i1, i2 = sc.calculate_interface_atoms(at1, at2, sc.INTERFACE_DISTANCE)
    combined = i1 + i2
    spacing = sc._choose_grid_spacing(combined, sc.PROBE_RADIUS, sc.GRID_SPACING)
    s1, n1 = sc._build_ses_surface(i1, sc.PROBE_RADIUS, spacing)
    s2, n2 = sc._build_ses_surface(i2, sc.PROBE_RADIUS, spacing)
    complex_surface, _ = sc._build_ses_surface(combined, sc.PROBE_RADIUS, spacing)
    thr = spacing * 1.25
    b1, bn1, a1, _ = sc._classify_surface_points(s1, n1, complex_surface, thr)
    b2, bn2, a2, _ = sc._classify_surface_points(s2, n2, complex_surface, thr)
    t1, _ = sc._trim_peripheral_band(b1, bn1, a1, sc.TRIM)
    t2, _ = sc._trim_peripheral_band(b2, bn2, a2, sc.TRIM)
    return {
        "interface_atoms_1": len(i1),
        "interface_atoms_2": len(i2),
        "spacing": spacing,
        "surface_1": len(s1),
        "surface_2": len(s2),
        "buried_1": len(b1),
        "buried_2": len(b2),
        "accessible_1": len(a1),
        "accessible_2": len(a2),
        "trimmed_1": len(t1),
        "trimmed_2": len(t2),
    }


df = pd.read_csv(CSV)
df["signed_error"] = df["python_sc"] - df["legacy_sc"]

under = df.sort_values("signed_error").head(5).copy()
over = df.sort_values("signed_error", ascending=False).head(5).copy()

rows = []
for label, subset in [("under", under), ("over", over)]:
    for _, row in subset.iterrows():
        stats = surface_stats(filename_to_model_path(row["filename"]))
        rows.append({
            "kind": label,
            "filename": row["filename"],
            "legacy_sc": row["legacy_sc"],
            "python_sc": row["python_sc"],
            "signed_error": row["signed_error"],
            **stats,
        })

out = pd.DataFrame(rows)
print(out.to_string(index=False))
