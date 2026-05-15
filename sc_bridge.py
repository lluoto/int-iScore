#!/usr/bin/env python3
"""SCASA bridge for sc_hts.

Single mode (C++ Accurate mode):
    sc_bridge.py <pdb_path> <chain1> <chain2>
    -> prints one numeric SC score to stdout for popen() compatibility.

Batch mode (C++ mode 2 Accurate optimization):
    sc_bridge.py --batch tasks.txt --output results.csv [--workers N]
    tasks.txt rows: pdb_path,chain1,chain2[,pdb_id]
    -> writes CSV rows in input order. Progress/errors go to stderr.
"""
from __future__ import annotations

# Prevent worker oversubscription before importing numpy/scipy stack.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover - tqdm may not be installed on minimal nodes
    def tqdm(iterable: Iterable, total: int | None = None, desc: str | None = None):
        return iterable

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "src"))

from int_iscore.utils.sc_calculator import calculate_sc_from_pdb

HEADER = ["PDB_ID", "Chains", "Mode", "TotalDots", "BuriedDots", "TrimmedDots", "Median_D", "Sc_Score"]


@dataclass(frozen=True)
class Task:
    index: int
    pdb_path: str
    chain1: str
    chain2: str
    pdb_id: str


def normalize_pid(pdb_path: str, pdb_id: str | None = None) -> str:
    if pdb_id:
        return Path(pdb_id.strip()).stem
    return Path(pdb_path).stem


def compute_task(task: Task) -> Tuple[int, List[object], str]:
    """Worker entry point. Returns (index, row, error). Never raises."""
    try:
        sc = calculate_sc_from_pdb(task.pdb_path, task.chain1, task.chain2)
        row = [task.pdb_id, f"{task.chain1}_{task.chain2}", "Accurate", 0, 0, 0, 0, float(sc)]
        return task.index, row, ""
    except Exception as exc:
        row = [task.pdb_id, f"{task.chain1}_{task.chain2}", "Accurate", 0, 0, 0, 0, -1.0]
        return task.index, row, f"{task.pdb_id}: {exc}"


def parse_tasks(path: str) -> List[Task]:
    tasks: List[Task] = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            if row[0].strip() == "pdb_path" or row[0].strip().startswith("PDB"):
                continue
            if len(row) < 3:
                print(f"[WARN] skip malformed task row: {row}", file=sys.stderr)
                continue
            pdb_path = row[0].strip()
            chain1 = row[1].strip()
            chain2 = row[2].strip()
            pdb_id = normalize_pid(pdb_path, row[3] if len(row) > 3 else None)
            tasks.append(Task(len(tasks), pdb_path, chain1, chain2, pdb_id))
    return tasks


def default_workers() -> int:
    cpus = os.cpu_count() or 2
    return max(1, cpus - 2)


def run_single(pdb_path: str, chain1: str, chain2: str) -> int:
    try:
        sc = calculate_sc_from_pdb(pdb_path, chain1, chain2)
        print(sc)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def run_batch(tasks_file: str, output_file: str, workers: int | None) -> int:
    tasks = parse_tasks(tasks_file)
    if not tasks:
        print("[WARN] no tasks found", file=sys.stderr)
        with open(output_file, "w", newline="") as handle:
            csv.writer(handle).writerow(HEADER)
        return 0

    n_workers = workers if workers and workers > 0 else default_workers()
    n_workers = min(n_workers, len(tasks))
    print(f"[INFO] Accurate batch: {len(tasks)} jobs, workers={n_workers}", file=sys.stderr)

    rows: List[List[object] | None] = [None] * len(tasks)
    errors = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        future_to_task = {pool.submit(compute_task, task): task for task in tasks}
        for future in tqdm(as_completed(future_to_task), total=len(future_to_task), desc="SCASA Accurate"):
            task = future_to_task[future]
            try:
                idx, row, err = future.result()
            except Exception as exc:  # defensive; compute_task should catch its own errors
                idx = task.index
                row = [task.pdb_id, f"{task.chain1}_{task.chain2}", "Accurate", 0, 0, 0, 0, -1.0]
                err = f"{task.pdb_id}: worker crashed: {exc}"
            rows[idx] = row
            if err:
                errors += 1
                print(f"[ERROR] {err}", file=sys.stderr)

    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        for row in rows:
            if row is not None:
                writer.writerow(row)

    print(f"[INFO] Accurate batch done: {len(tasks) - errors} ok, {errors} failed", file=sys.stderr)
    return 0 if errors == 0 else 2


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="SCASA bridge for sc_hts")
    parser.add_argument("pdb_path", nargs="?", help="PDB file for single-job mode")
    parser.add_argument("chain1", nargs="?", help="chain 1 for single-job mode")
    parser.add_argument("chain2", nargs="?", help="chain 2 for single-job mode")
    parser.add_argument("--batch", help="batch task file: pdb_path,chain1,chain2[,pdb_id]")
    parser.add_argument("--output", default="results.csv", help="batch output CSV path")
    parser.add_argument("--workers", type=int, default=0, help="worker count (default: os.cpu_count()-2)")
    args = parser.parse_args(argv)

    if args.batch:
        return run_batch(args.batch, args.output, args.workers)
    if args.pdb_path and args.chain1 and args.chain2:
        return run_single(args.pdb_path, args.chain1, args.chain2)
    parser.print_usage(sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())