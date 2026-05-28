#!/usr/bin/env bash
set -euo pipefail

# High-throughput SC pipeline: Accurate batch scoring.
# Input CSV format: pdb_path,chain1,chain2[,pdb_id]
# All jobs scored with AccurateConnolly (SCASA bridge, CCP4-grade, +/-2% accuracy).

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input_jobs.csv> <output_dir> [workers=auto]" >&2
  echo "  input_jobs.csv rows: pdb_path,chain1,chain2[,pdb_id]" >&2
  echo "  output files: accurate_results.csv" >&2
  exit 1
fi

INPUT_CSV="$1"
OUT_DIR="$2"
WORKERS="${3:-0}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SC_HTS="$SCRIPT_DIR/sc_hts"

mkdir -p "$OUT_DIR"
ACCURATE_JOBS="$OUT_DIR/accurate_jobs.csv"
ACCURATE_RESULTS="$OUT_DIR/accurate_results.csv"
ACCURATE_LOG="$OUT_DIR/accurate.log"

if [[ ! -x "$SC_HTS" ]]; then
  echo "[INFO] Building sc_hts..." >&2
  g++ -std=c++17 -O3 -fopenmp -I"$SCRIPT_DIR" "$SCRIPT_DIR/sc_hts.cpp" -o "$SC_HTS"
fi

echo "[STEP] Preparing Accurate batch jobs..." >&2
python3 - "$INPUT_CSV" "$ACCURATE_JOBS" <<'PY'
import csv, sys, pathlib
inp, out = sys.argv[1], sys.argv[2]
with open(inp, newline='') as fin, open(out, 'w', newline='') as fout:
    writer = csv.writer(fout)
    writer.writerow(['pdb_path','chain1','chain2','mode','pdb_id'])
    for row in csv.reader(fin):
        if not row or row[0].startswith('#') or row[0] in ('pdb_path','PDB_ID'):
            continue
        pdb = row[0].strip()
        c1 = row[1].strip()
        c2 = row[2].strip()
        pid = row[3].strip() if len(row) > 3 and row[3].strip() else pathlib.Path(pdb).stem
        writer.writerow([pdb, c1, c2, 1, pid])
PY

echo "[STEP] AccurateConnolly batch scoring ($ACCURATE_JOBS -> $ACCURATE_RESULTS)" >&2
"$SC_HTS" "$ACCURATE_JOBS" ? ? 2 > "$ACCURATE_RESULTS" 2> "$ACCURATE_LOG"

echo "Done. Results: $ACCURATE_RESULTS" >&2
echo "Log: $ACCURATE_LOG" >&2
