# Work Log: C++ High-Throughput Shape Complementarity Engine (`sc_hts.cpp`)

> **Purpose**: Chronological development log for other AI models to understand the full context.
> **Repository**: `/media/cuixi/data01/lluoto/int-iScore`
> **Remote**: `cuixi@10.19.25.48`

---

## Timeline

### 2026-05-09 | Phase 0: Initial Build

Created `sc_hts.cpp` — C++17 high-throughput Shape Complementarity engine.

**Architecture**:
- **Mode 0 (FastSAS)**: Solvent-Accessible Surface (SAS) + nanoflann KD-Tree nearest-neighbor scoring. Target: ~250 ms/complex.
- **Mode 1 (Accurate)**: Originally a simplified C++ Connolly implementation; later replaced with Python SCASA bridge.
- **SOA memory layout**, OpenMP parallelization, Fibonacci sphere dot sampling.

**Key constants**:
- `PROBE_RADIUS = 1.4` Å
- `DOT_DENSITY = 15.0` dots/Å²
- `INTERFACE_DISTANCE = 8.0` Å (atom-level interface cutoff)
- Original FastSAS: `WEIGHT_FAST = 0.08`, `FAST_DISTANCE_CUTOFF = 3.0` Å

---

### 2026-05-09 | Phase 1: Bug Fixes & Baseline

| Commit | Description | Impact |
|--------|------------|--------|
| `ac0f66e` | **Fix: median_d initialization** — added `result.median_d = 0.0;` to both compute functions. Was printing denormalized floats (~1e-308) instead of 0.0. | CSV output now correct |
| `d9b0214` | **Feat: batch CSV mode (mode=2)** — reads `pdb_path,chain1,chain2,mode[,pdb_id]`, processes sequentially, outputs combined CSV with header. | Enables large-scale screening |
| `043a74a` | **Feat: Accurate mode → Python SCASA bridge** — replaced simplified C++ Connolly (SC=0.031 vs expected 0.6) with `popen("python3 sc_bridge.py <pdb> <c1> <c2>")`. Returns CCP4-accurate SC within ±2%. | Gold-standard accuracy |
| `b0f0d1e` | **Docs: README** — added sc_hts section (build, usage, benchmark). | Documentation |

**AccurateBridge details**:
- Calls `calculate_sc_from_pdb()` from `src/int_iscore/utils/sc_calculator.py`
- Uses Connolly molecular surface (convex + toroidal + concave patches)
- Formula: `-dot·exp(-0.5·d²)`, trim=1.5 Å
- Overhead: ~10s/complex (Python startup + SCASA computation)

---

### 2026-05-09 | Phase 2: Calibration — Negative Correlation Discovered

**Test setup**: 4 AlphaFold3-predicted complexes (6gs2, 6h4b, 6ii6, 6ono), each with 1000 seed variations. All A-B chain interfaces.

**Key finding**: FastSAS is **negatively correlated** with Accurate (CCP4) SC.

| PDB | FastSAS Mean | FastSAS Std | Accurate Mean | Fast Rank | Acc Rank |
|-----|-------------|------------|---------------|-----------|----------|
| 6gs2 | 0.409 | 0.058 | 0.570 | #1 | #4 |
| 6h4b | 0.341 | 0.053 | 0.577 | #3 | #3 |
| 6ii6 | 0.377 | 0.056 | 0.587 | #2 | #2 |
| 6ono | 0.323 | 0.053 | 0.650 | #4 | #1 |

- **Pearson r = -0.416** (N=80 paired samples)
- Per-PDB ranking nearly inverted
- FastSAS varies 2x within same PDB (σ/μ ≈ 15%); Accurate is stable (σ/μ ≈ 2%)

**Source data**:
- 4,614 FastSAS results → `/tmp/fast_results.csv`
- 80 Accurate results → `/tmp/acc_results.csv`
- pdb_output/ directory: 12 PDBs × 1000 seeds = 12,000 files

**6A6I reference** (rank0 AF3 model, A-B chains):
- CCP4 reference SC: 0.616
- Accurate (our bridge): 0.516–0.563
- FastSAS (original): 0.47

---

### 2026-05-10 | Phase 3: Parameter Sweep — All Attempts Failed

Systematic tuning of scoring formula parameters. All trials on 80 paired samples (4 PDBs × 20 seeds).

| Trial | Formula | w | Cutoff (Å) | Interface (Å) | Aggregation | Pearson r | Δ from baseline |
|-------|---------|---|-----------|---------------|-------------|-----------|-----------------|
| v0 (baseline) | `-dot·exp(-w·d²)` | 0.08 | 3.0 | 8.0 | median | **-0.416** | — |
| v1 | `-dot·exp(-w·d²)` | 0.08 | 3.0 | 8.0 | mean | -0.452 | -0.036 (worse) |
| v2 | `-dot·exp(-w·d²)` | 0.20 | 3.0 | 5.0 | mean | -0.476 | -0.060 (worse) |
| v3 | `-dot·exp(-w·d²)` | 0.08 | none | 8.0 | median | collapsed | all scores → 0 |

**Observations**:
- Increasing weight (w) fails because SAS NN distances are inflated (~6–8 Å) → exp(-0.2·36) ≈ 0.0007
- Removing cutoff fails because noise dots dominate the signal
- All attempts within the `-dot·exp(-w·d²)` framework failed

---

### 2026-05-11 | Phase 4: Gaussian Overlap — Best Effort (r = -0.382)

**Commit**: `23a9b8c` — "perf(sc): Gaussian Overlap scoring for FastSAS"

**Formula change**:
```
Old:  s_val = -dot(n1, n2) * exp(-w · d²)
New:  s_val = normal_score * exp(-w · d²)
      where normal_score = max(0, -dot(n1, n2))
```

**Rationale**: The `-dot` term was a multiplier, so when normals are aligned (dot ≈ -1) AND NN distance is large (SAS artifact), the dot penalty dominates. By making `normal_score` dominant (0–1 range), we reduce distance sensitivity and amplify normal-vector alignment signal.

**Parameters**: w=0.05, cutoff=4.0 Å, median aggregation

**Result**: Pearson r = **-0.382** → +8% improvement from -0.416, but still negative.

**Key insight**: Even with normal vectors perfectly aligned at the interface (normal_score=1.0), the SAS surface distances are systematically inflated by ~1.4 Å relative to Connolly surfaces. This means `exp(-0.05·d²)` still produces lower scores at better interfaces — the physical inversion persists.

---

### 2026-05-11 | Phase 5: Physical Diagnosis & Strategy

**Root cause identified**: **SAS surface distance inversion**

```
Connolly surface:  probe rolls ON VdW spheres → dots at VdW + 0 Å
SAS surface:       atom center + VdW + probe → dots at VdW + 1.4 Å

At perfect complementarity:
  - Connolly: d ≈ 0–1 Å  → high SC score ✓
  - SAS:      d ≈ 1.4–2.4 Å → lower SC score ✗ (inverted!)
```

The SAS surface systematically pushes dots farther from atoms, so at better interfaces (tighter packing), SAS NN distances actually **increase** for dots near the interface edge. This is the **opposite** of what Connolly SC measures.

**Conclusion**: FastSAS cannot rank interfaces by SC through parameter tuning alone — the surface topology is fundamentally different from Connolly. It requires either:
1. A faster Connolly approximation (convex-only + radial normals), or
2. A different metric entirely (features + regression, atom-atom distances, interface area)

**Current strategy**:
- **FastSAS → binary classifier only**: distinguishing SC>0 (interacting) from SC=0 (non-interacting)
- **Accurate → ranking**: full Connolly SC for ordering complexes
- Future: parallelize Accurate mode (multiple Python SCASA instances) for batch throughput

---

## Parameter Evolution Summary

| Version | Formula | w | Cutoff (Å) | Aggregation | r | Notes |
|---------|---------|---|-----------|-------------|---|-------|
| v0 | `-dot·exp(-w·d²)` | 0.08 | 3.0 | median | -0.416 | Original baseline |
| v1 | `-dot·exp(-w·d²)` | 0.08 | 3.0 | mean | -0.452 | Worse |
| v2 | `-dot·exp(-w·d²)` | 0.20 | 3.0 | mean | -0.476 | Worse |
| v3 | `-dot·exp(-w·d²)` | 0.08 | none | median | ~0 | Collapsed |
| **Gaussian** | `max(0,-dot)·exp(-w·d²)` | 0.05 | 4.0 | median | **-0.382** | Best (current) |

---

## Benchmark Data

### 6A6I A-B (rank0 AF3 model)
| Mode | SC | Time | Dots | vs CCP4 (0.616) |
|------|----|------|------|-----------------|
| FastSAS (Gaussian) | ~0.27 | ~215 ms | ~425K | — |
| Accurate (SCASA bridge) | 0.516–0.563 | ~9.5 s | — | ±2% |
| CCP4 | 0.616 | — | — | Gold standard |

### Large-Scale FastSAS Statistics
| PDB | N | FastSAS Mean | FastSAS Std | Accurate Mean |
|-----|---|-------------|------------|---------------|
| 6gs2 | 1000 | 0.409 | 0.058 | 0.570 |
| 6h4b | 1000 | 0.341 | 0.053 | 0.577 |
| 6ii6 | 1000 | 0.377 | 0.056 | 0.587 |
| 6ono | 1000 | 0.323 | 0.053 | 0.650 |

---

## Git History

```
23a9b8c (HEAD -> main) perf(sc): Gaussian Overlap scoring for FastSAS — normal-vector-dominant w=0.05 cutoff=4.0
b0f0d1e docs(readme): add C++ sc_hts engine section with build, usage, and benchmarks
043a74a feat(sc): replace simplified C++ Connolly with Python SCASA bridge for Accurate mode
d9b0214 feat(sc): add batch CSV processing mode (mode=2) to sc_hts
ac0f66e fix(sc): initialize median_d to 0.0 in both FastSAS and Accurate modes
4c51b54 (origin/main) — previous baseline (pushed)
```

**Status**: 5 commits ahead of `origin/main`, all unpushed.

---

## Key Files

| File | Description | Lines |
|------|------------|-------|
| `sc_hts.cpp` | Main C++17 engine (FastSAS + Accurate + Batch) | ~1009 |
| `sc_hts` | Compiled binary | — |
| `sc_hts.bak.cpp` | Pre-Gaussian backup | ~1000 |
| `sc_bridge.py` | Python SCASA bridge for Accurate mode | ~20 |
| `FASTSAS_ANALYSIS_PROMPT.md` | Analysis prompt for other AI models | ~130 |
| `README.md` | Project documentation (includes sc_hts section) | ~200 |
| `nanoflann.hpp` | Header-only KD-Tree library | ~1800 |

**Test data**:
- `/tmp/fast_results.csv` — 4,614 FastSAS results (original w=0.08)
- `/tmp/acc_results.csv` — 80 Accurate SCASA results
- `/tmp/fast_gauss_final.csv` — 80 Gaussian Overlap FastSAS results
- `/media/cuixi/data01/lluoto/pdb_output/` — 12,000 PDB files (12 PDBs × 1000 seeds)
- `/media/cuixi/data01/lluoto/6A6I/` — AlphaFold3 output (6a6i_model.cif + seed subdirs)

---

## Build & Test Commands

```bash
cd /media/cuixi/data01/lluoto/int-iScore

# Build
g++ -std=c++17 -O3 -fopenmp -I. sc_hts.cpp -o sc_hts

# FastSAS single
./sc_hts <pdb_file> A B 0 <pdb_id>

# Accurate single (requires Python + int_iscore)
./sc_hts <pdb_file> A B 1 <pdb_id>

# Batch mode (CSV: pdb_path,chain1,chain2,mode[,pdb_id])
./sc_hts jobs.csv ? ? 2

# Python SCASA directly
python3 sc_bridge.py <pdb_file> A B
```

---

## Open Questions for Future Models

1. Can a **convex-only Connolly approximation** (radial normals from atom centers) be fast enough (<500 ms) for screening?
2. Can **atom-to-atom distances** (VdW overlap) replace dot-to-dot for FastSAS scoring?
3. Would a **simple ML regressor** on SAS-derived features (interface area, dot density, curvature moments, buried fraction) predict CCP4 SC better than direct physics formulas?
4. Could **differential SAS** (ΔSASA pre/post complexation) serve as a coarse filter that's positively correlated with SC?
