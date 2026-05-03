import math
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, "/media/cuixi/data01/lluoto/int-iScore/src")
sys.path.insert(0, "/media/cuixi/data01/lluoto/int-iScore")

from int_iscore.utils.sc_calculator import calculate_sc_from_pdb


BASE = Path("/media/cuixi/data01/lluoto/6A6I")
LEGACY = Path("/media/cuixi/data01/lluoto/int-iScore/all_3_6.csv")
OUT = Path("/media/cuixi/data01/lluoto/int-iScore/int_iscore_temp/6A6I_sc_comparison.csv")


def filename_to_model_path(name: str) -> Path:
    suffix = name[len("6A6I_"):]
    if not suffix.endswith("model"):
        raise ValueError(name)
    stem = suffix[:-len("model")]
    return BASE / stem / "model.pdb"


def main() -> None:
    df = pd.read_csv(LEGACY)
    rows = df[df["filename"].astype(str).str.startswith("6A6I_")].copy()
    results = []
    for _, row in rows.iterrows():
        filename = row["filename"]
        model_path = filename_to_model_path(filename)
        exists = model_path.exists()
        py_sc = math.nan
        error = ""
        if exists:
            try:
                py_sc = calculate_sc_from_pdb(str(model_path), "A", "B")
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
        else:
            error = "missing model.pdb"
        legacy_sc = float(row["sc"])
        results.append(
            {
                "filename": filename,
                "model_path": str(model_path),
                "exists": exists,
                "legacy_sc": legacy_sc,
                "python_sc": py_sc,
                "abs_error": abs(py_sc - legacy_sc) if exists and not math.isnan(py_sc) else math.nan,
                "rel_error": abs(py_sc - legacy_sc) / legacy_sc if exists and not math.isnan(py_sc) and legacy_sc != 0 else math.nan,
                "error": error,
            }
        )
    out_df = pd.DataFrame(results)
    out_df.to_csv(OUT, index=False)

    valid = out_df[out_df["python_sc"].notna()].copy()
    print(f"rows={len(out_df)} valid={len(valid)} missing_or_failed={len(out_df)-len(valid)}")
    if len(valid):
        print(f"mean_legacy={valid['legacy_sc'].mean():.6f}")
        print(f"mean_python={valid['python_sc'].mean():.6f}")
        print(f"mae={valid['abs_error'].mean():.6f}")
        print(f"median_abs_error={valid['abs_error'].median():.6f}")
        print(f"rmse={(valid['abs_error'].pow(2).mean() ** 0.5):.6f}")
        print(f"mean_rel_error={valid['rel_error'].mean():.6f}")
        within_005 = (valid['abs_error'] <= 0.05).mean()
        within_01 = (valid['abs_error'] <= 0.10).mean()
        print(f"within_0.05={within_005:.6f}")
        print(f"within_0.10={within_01:.6f}")
        print("worst5")
        print(valid.sort_values('abs_error', ascending=False)[['filename','legacy_sc','python_sc','abs_error','rel_error']].head(5).to_string(index=False))
    print(f"wrote={OUT}")


if __name__ == "__main__":
    main()
