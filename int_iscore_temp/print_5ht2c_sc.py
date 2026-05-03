import pandas as pd
import sys

sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore/src")
sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore")

from int_iscore.utils import calculate_sc


legacy = pd.read_csv("/mnt/sdb2/lluoto/int-iScore/5HT2C_intiscore.csv").head(5)

for _, row in legacy.iterrows():
    name = row["filename"]
    pdb = f"/mnt/sdb2/lluoto/5-HT2C-2CD28/{name.replace('.cif', '.pdb')}"
    current = calculate_sc(pdb, ["A"], ["B", "C"])["sc_mean"]
    legacy_sc = row["sc"] if "sc" in row else "NA"
    print(f"{name}\tlegacy_sc={legacy_sc}\tcurrent_sc={current}")
