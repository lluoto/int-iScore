import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore/src")
sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore")

from int_iscore.utils import calculate_sc
from int_iscore.utils.sc_calculator import calculate_sc_from_pdb


ROOT = Path("/mnt/sdb2/lluoto/5-HT2C-2CD28")
HELPER = Path("/mnt/sdb2/lluoto/int-iScore/int_iscore_temp/run_ccp4_sc.sh")


def run_ccp4(pdb: Path, *chains: str) -> float:
    proc = subprocess.run(
        ["bash", str(HELPER), str(pdb), *chains],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in proc.stdout.splitlines():
        if "Shape complementarity statistic Sc =" in line:
            return float(line.split("=")[-1].strip())
    raise RuntimeError(f"SC line not found for {pdb} {chains}\n{proc.stdout}\n{proc.stderr}")


for i in range(5):
    pdb = ROOT / f"pred.rank_{i}.pdb"
    py_ab = calculate_sc_from_pdb(str(pdb), "A", "B")
    py_ac = calculate_sc_from_pdb(str(pdb), "A", "C")
    py_abc = calculate_sc(str(pdb), ["A"], ["B", "C"])["sc_mean"]

    ccp4_ab = run_ccp4(pdb, "A", "B")
    ccp4_ac = run_ccp4(pdb, "A", "C")
    ccp4_abc = run_ccp4(pdb, "A", "B", "C")

    print(f"pred.rank_{i}")
    print(f"  A-B     CCP4={ccp4_ab:.6f}  PY={py_ab:.6f}")
    print(f"  A-C     CCP4={ccp4_ac:.6f}  PY={py_ac:.6f}")
    print(f"  A-(B+C) CCP4={ccp4_abc:.6f}  PY={py_abc:.6f}")
