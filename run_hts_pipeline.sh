#!/usr/bin/env bash
set -euo pipefail

# High-throughput SC pipeline: FastSAS_DEPRECATED prefilter + Accurate parallel scoring.
# Input CSV format: pdb_path,chain1,chain2[,mode,pdb_id]
# FastSAS is DEPRECATED for cross-PDB ranking; use only within the same target/PDB seed set.

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input_jobs.csv> <output_dir> [top_fraction=0.5] [workers=auto]" >&2
  echo "  input_jobs.csv rows: pdb_path,chain1,chain2[,mode,pdb_id]" >&2
  echo "  output files: fastsas_results.csv, top_jobs.csv, accurate_results.csv" >&2
  exit 1
fi

INPUT_CSV="$1"
OUT_DIR="$2"
TOP_FRAC="${3:-0.5}"
WORKERS="${4:-0}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SC_HTS="$SCRIPT_DIR/sc_hts"

mkdir -p "$OUT_DIR"
FAST_JOBS="$OUT_DIR/fastsas_jobs.csv"
FAST_RESULTS="$OUT_DIR/fastsas_results.csv"
FAST_LOG="$OUT_DIR/fastsas.log"
TOP_JOBS="$OUT_DIR/top_jobs.csv"
ACCURATE_RESULTS="$OUT_DIR/accurate_results.csv"
ACCURATE_LOG="$OUT_DIR/accurate.log"

if [[ ! -x "$SC_HTS" ]]; then
  echo "[INFO] Building sc_hts..." >&2
  g++ -std=c++17 -O3 -fopenmp -I"$SCRIPT_DIR" "$SCRIPT_DIR/sc_hts.cpp" -o "$SC_HTS"
fi

echo "[STEP 1] FastSAS_DEPRECATED coarse filter" >&2
python3 - "$INPUT_CSV" "$FAST_JOBS" <<'PY'
import csv, sys
inp, out = sys.argv[1], sys.argv[2]
with open(inp, newline='') as fin, open(out, 'w', newline='') as fout:
    reader = csv.reader(fin)
    writer = csv.writer(fout)
    writer.writerow(['pdb_path','chain1','chain2','mode','pdb_id'])
    for row in reader:
        if not row or row[0].startswith('#') or row[0] in ('pdb_path','PDB_ID'):
            continue
        pdb = row[0].strip(); c1 = row[1].strip(); c2 = row[2].strip()
        pid = row[4].strip() if len(row) > 4 and row[4].strip() else pdb
        writer.writerow([pdb, c1, c2, 0, pid])
PY
"$SC_HTS" "$FAST_JOBS" ? ? 2 > "$FAST_RESULTS" 2> "$FAST_LOG"

echo "[STEP 2] Selecting top fraction: $TOP_FRAC" >&2
python3 - "$INPUT_CSV" "$FAST_RESULTS" "$TOP_JOBS" "$TOP_FRAC" <<'PY'
import csv, sys, pathlib
input_csv, fast_csv, top_csv, frac_s = sys.argv[1:5]
frac = float(frac_s)
orig = {}
with open(input_csv, newline='') as f:
    for row in csv.reader(f):
        if not row or row[0].startswith('#') or row[0] in ('pdb_path','PDB_ID'):
            continue
        pdb, c1, c2 = row[0].strip(), row[1].strip(), row[2].strip()
        pid = row[4].strip() if len(row) > 4 and row[4].strip() else pathlib.Path(pdb).stem
        orig[pathlib.Path(pid).stem] = (pdb, c1, c2, pathlib.Path(pid).stem)

scores = []
with open(fast_csv, newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        if not row or row[0] == 'PDB_ID':
            continue
        try:
            pid = pathlib.Path(row[0].strip()).stem
            sc = float(row[-1])
            if pid in orig:
                scores.append((sc, pid))
        except Exception:
            pass
scores.sort(reverse=True)
keep = max(1, int(len(scores) * frac + 0.999999))
with open(top_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['pdb_path','chain1','chain2','mode','pdb_id'])
    for sc, pid in scores[:keep]:
        pdb, c1, c2, out_pid = orig[pid]
        writer.writerow([pdb, c1, c2, 1, out_pid])
print(f'[INFO] selected {keep}/{len(scores)} jobs', file=sys.stderr)
PY

echo "[STEP 3] Accurate parallel batch" >&2
if [[ "$WORKERS" == "0" || "$WORKERS" == "auto" ]]; then
  "$SC_HTS" "$TOP_JOBS" ? ? 2 > "$ACCURATE_RESULTS" 2> "$ACCURATE_LOG"
else
  # Direct bridge worker override for advanced users.
  TASKS="$OUT_DIR/accurate_tasks.csv"
  python3 - "$TOP_JOBS" "$TASKS" <<'PY'
import csv, sys
inp, out = sys.argv[1], sys.argv[2]
with open(inp, newline='') as fin, open(out, 'w', newline='') as fout:
    for row in csv.reader(fin):
        if not row or row[0] in ('pdb_path','PDB_ID') or row[0].startswith('#'):
            continue
        fout.write(','.join([row[0], row[1], row[2], row[4] if len(row) > 4 else row[0]]) + '\n')
PY
  python3 "$SCRIPT_DIR/sc_bridge.py" --batch "$TASKS" --output "$ACCURATE_RESULTS" --workers "$WORKERS" 2> "$ACCURATE_LOG"
fi

echo "[DONE]" >&2
echo "  FastSAS:  $FAST_RESULTS" >&2
echo "  Top jobs: $TOP_JOBS" >&2
echo "  Accurate: $ACCURATE_RESULTS" >&2