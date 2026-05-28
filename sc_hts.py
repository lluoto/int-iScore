#!/usr/bin/env python3
"""
sc_hts.py — 纯 Python 高通量 Shape Complementarity 批次调度引擎

替代原 C++ 前端 (sc_hts.cpp)，消除跨语言 popen 调用和中间临时文件。
直接读取 CSV → 内存路由分发 → multiprocessing(spawn) 并行计算 → 结果汇总输出。

用法:
    # 单 PDB 模式（兼容原 C++ popen 调用）
    python3 sc_hts.py <pdb_file> <chain1> <chain2>

    # 批次模式
    python3 sc_hts.py --batch input.csv --output results.csv [--workers N]

CSV 格式:
    pdb_path,chain1,chain2,mode[,pdb_id]
    mode: 0/1=AccurateConnolly（兼容列，实际都走 Accurate）

环境要求: Python 3.10+, int_iscore 包已安装

设计要点:
    - spawn 上下文: int_iscore 底层 C++ 扩展非 fork-safe，必须用 spawn
    - maxtasksperchild=50: 防止 C++ 扩展内存泄漏累积，定期回收 Worker
    - 容错: 单个 PDB 崩溃不影响整体批次
    - 零中间文件: 输入 CSV 直接解析，输出直接写入，无临时文件
"""

from __future__ import annotations

import os
import sys

# ── 环境锁（必须在任何 numpy/scipy 导入之前） ──────────────────────────
# 防止底层 C++ 扩展的线程库在线程池中相互抢占
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import csv
import multiprocessing
import signal
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ── 路径设置：确保可以 import int_iscore ───────────────────────────────
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from int_iscore.utils.sc_calculator import calculate_sc_from_pdb

# CSV 输出表头（兼容原 C++ 格式）
HEADER = ["PDB_ID", "Chains", "Mode", "TotalDots", "BuriedDots", "TrimmedDots", "Median_D", "Sc_Score"]


# ══════════════════════════════════════════════════════════════════════════════
# 数据结构
# ══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Task:
    """单个计算任务。index 用于保持输出顺序与输入一致。"""
    index: int
    pdb_path: str
    chain1: str
    chain2: str
    pdb_id: str


# ══════════════════════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════════════════════

def normalize_pid(pdb_path: str, pdb_id: Optional[str] = None) -> str:
    """从路径提取 PDB ID（去掉目录和扩展名）。"""
    if pdb_id and pdb_id.strip():
        return Path(pdb_id.strip()).stem
    return Path(pdb_path).stem


def parse_tasks(csv_path: str) -> List[Task]:
    """
    解析输入 CSV，返回 Task 列表。
    CSV 列: pdb_path, chain1, chain2, mode[, pdb_id]
    跳过表头行和空行，# 开头的注释行。
    """
    tasks: List[Task] = []
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            # 跳过空行、注释行
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            # 跳过表头和 PDB_ID 开头的行（兼容结果 CSV 被误传）
            col0 = row[0].strip()
            if col0 in ("pdb_path", "PDB_ID") or col0.startswith("PDB"):
                continue
            if len(row) < 3:
                print(f"[WARN] skip malformed row (need >=3 cols): {row}", file=sys.stderr)
                continue

            pdb_path = row[0].strip()
            chain1 = row[1].strip()
            chain2 = row[2].strip()
            # mode 列兼容：有就保留，无则默认 Accurate
            # pdb_id 列兼容：有则用，无则从路径提取
            pdb_id = normalize_pid(pdb_path, row[4] if len(row) > 4 else row[3] if len(row) > 3 else None)

            tasks.append(Task(
                index=len(tasks),
                pdb_path=pdb_path,
                chain1=chain1,
                chain2=chain2,
                pdb_id=pdb_id,
            ))
    return tasks


def compute_task(task: Task) -> Tuple[int, List[Any], str]:
    """
    Worker 入口 — 单个 PDB 的 SC 计算。

    必须定义为模块级函数（multiprocessing spawn 要求可 picklable）。

    Returns:
        (index, row, error_message)
        - index: 任务序号，用于结果排序
        - row: 输出 CSV 行列表（8 列）
        - error_message: 空字符串表示成功，否则为错误描述
    """
    try:
        # ── 核心调用：int_iscore C++ 扩展 ──
        sc_score = calculate_sc_from_pdb(task.pdb_path, task.chain1, task.chain2)
        row = [
            task.pdb_id,
            f"{task.chain1}_{task.chain2}",
            "Accurate",
            0, 0, 0, 0,           # Total/Buried/Trimmed/Median — Accurate 模式不计算这些
            float(sc_score),
        ]
        return task.index, row, ""
    except FileNotFoundError:
        row = [task.pdb_id, f"{task.chain1}_{task.chain2}", "Accurate", 0, 0, 0, 0, -1.0]
        return task.index, row, f"{task.pdb_id}: PDB file not found: {task.pdb_path}"
    except ValueError as exc:
        row = [task.pdb_id, f"{task.chain1}_{task.chain2}", "Accurate", 0, 0, 0, 0, -1.0]
        return task.index, row, f"{task.pdb_id}: chain or atom error: {exc}"
    except Exception as exc:
        row = [task.pdb_id, f"{task.chain1}_{task.chain2}", "Accurate", 0, 0, 0, 0, -1.0]
        return task.index, row, f"{task.pdb_id}: {type(exc).__name__}: {exc}"


# ══════════════════════════════════════════════════════════════════════════════
# 初始化 Worker（spawn 上下文）
# ══════════════════════════════════════════════════════════════════════════════

def _worker_init():
    """Worker 进程初始化 — 忽略 SIGINT，由主进程统一管理。"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


# ══════════════════════════════════════════════════════════════════════════════
# 核心调度
# ══════════════════════════════════════════════════════════════════════════════

def default_workers() -> int:
    """默认 Worker 数 = CPU - 2（最少 1）。"""
    cpus = os.cpu_count() or 2
    return max(1, cpus - 2)


def run_batch(input_csv: str, output_csv: str, workers: int = 0) -> int:
    """
    批次模式主入口。

    流程:
        1. 解析 CSV → Task 列表（内存中）
        2. 创建 spawn Pool + maxtasksperchild=50
        3. map_async 分发所有任务
        4. 收集结果，按 index 排序
        5. 写入输出 CSV

    Returns:
        0: 全部成功
        2: 有部分失败
    """
    # ── Step 1: 解析输入 CSV ──
    tasks = parse_tasks(input_csv)
    if not tasks:
        print("[WARN] no valid tasks found in input CSV", file=sys.stderr)
        # 写空表头
        with open(output_csv, "w", newline="", encoding="utf-8") as handle:
            csv.writer(handle).writerow(HEADER)
        return 0

    n_workers = workers if workers > 0 else default_workers()
    n_workers = min(n_workers, len(tasks))
    print(f"[INFO] Batch: {len(tasks)} jobs, workers={n_workers}, maxtasksperchild=50", file=sys.stderr)

    # ── Step 2: 创建 spawn Pool ──
    # spawn 上下文: 确保 C++ 扩展在干净进程中初始化，避免 fork 后遗留的锁和状态
    # maxtasksperchild=50: 定期回收 Worker，防止 C++ 扩展微内存泄漏累积
    ctx = multiprocessing.get_context("spawn")
    t_start = time.monotonic()

    results_map: Dict[int, List[Any]] = {}
    errors: List[str] = []

    try:
        with ctx.Pool(
            processes=n_workers,
            initializer=_worker_init,
            maxtasksperchild=50,
        ) as pool:
            # 使用 starmap_async 分发任务
            async_result = pool.starmap_async(
                compute_task,
                [(task,) for task in tasks],
                chunksize=max(1, len(tasks) // (n_workers * 4)),  # 动态分块减少 IPC
            )

            # 等待完成，带超时（None = 无限等待）
            # 同时允许 Ctrl+C 中断
            while not async_result.ready():
                async_result.wait(timeout=1.0)

            # 收集结果
            for idx, row, err in async_result.get():
                results_map[idx] = row
                if err:
                    errors.append(err)

    except KeyboardInterrupt:
        print("\n[ABORT] Interrupted by user", file=sys.stderr)
        pool.terminate()
        pool.join()
        return 130
    except Exception as exc:
        print(f"[FATAL] Pool execution failed: {exc}", file=sys.stderr)
        return 1

    elapsed = time.monotonic() - t_start

    # 打印错误
    for err in errors:
        print(f"[ERROR] {err}", file=sys.stderr)

    # ── Step 3: 写入输出 CSV ──
    # 按 index 排序以保证输出顺序与输入一致
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        for i in range(len(tasks)):
            if i in results_map:
                writer.writerow(results_map[i])
            else:
                # 兜底：任务未返回（理论上不会走到这里）
                t = tasks[i]
                writer.writerow([t.pdb_id, f"{t.chain1}_{t.chain2}", "Accurate", 0, 0, 0, 0, -1.0])

    n_ok = len(tasks) - len(errors)
    rate = len(tasks) / elapsed if elapsed > 0 else 0
    print(f"[INFO] Done: {n_ok}/{len(tasks)} ok, {len(errors)} failed, "
          f"{elapsed:.1f}s ({rate:.2f} jobs/s)", file=sys.stderr)

    return 0 if not errors else 2


def run_single(pdb_path: str, chain1: str, chain2: str) -> int:
    """单 PDB 模式 — 兼容原 C++ popen 调用。输出纯数字到 stdout。"""
    try:
        sc = calculate_sc_from_pdb(pdb_path, chain1, chain2)
        print(sc)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Pure-Python high-throughput Shape Complementarity batch engine",
    )
    # 单 PDB 模式（位置参数）
    parser.add_argument("pdb_path", nargs="?", help="PDB file path (single-job mode)")
    parser.add_argument("chain1", nargs="?", help="Chain 1 ID (single-job mode)")
    parser.add_argument("chain2", nargs="?", help="Chain 2 ID (single-job mode)")
    # 批次模式（命名参数）
    parser.add_argument("--batch", help="Input CSV: pdb_path,chain1,chain2,mode[,pdb_id]")
    parser.add_argument("--output", default="results.csv", help="Output CSV path (default: results.csv)")
    parser.add_argument("--workers", type=int, default=0,
                        help=f"Worker count (default: cpu_count-2 = {default_workers()})")

    args = parser.parse_args(argv)

    if args.batch:
        return run_batch(args.batch, args.output, args.workers)

    if args.pdb_path and args.chain1 and args.chain2:
        return run_single(args.pdb_path, args.chain1, args.chain2)

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
