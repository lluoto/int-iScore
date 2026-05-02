# Current Status: SC and BSA

## Worktree

- Active project path on server: `/media/cuixi/data01/lluoto/int-iScore`

## What was done so far

### BSA / freesasa

1. Verified that the **active package path** already uses `freesasa` for BSA in:
   - `src/int_iscore/core/metrics.py`
2. Verified that the old **VMD/Tcl BSA path** still existed in legacy scripts, especially:
   - `calculate_md.py`
3. Patched legacy scripts so they use the same selected-chain `freesasa` BSA idea as the active package:
   - `calculate_md.py`
   - `calculate_non_ref_af3.py`
   - `calculate_af3.py`

### SC / shape complementarity

1. Replaced the earlier simple atom-sphere scoring with a more CCP4-like Python path.
2. Added / refined:
   - CCP4-like defaults
   - CCP4-style atom radii table
   - combined-molecule handling for multi-chain partner SC
   - CCP4 comparison helper script:
     - `int_iscore_temp/run_ccp4_sc.sh`
3. Benchmarked against the original CCP4 `sc` program on representative cases.
4. Benchmarked repeatedly against `6A6I` and `5HT2C` datasets.

### Frustratometer / frustration score

1. Inspected legacy script:
   - `frustratometer2-master/run.py`
2. Confirmed it is a standalone batch script tied to old file paths and old CSV workflow, for example:
   - reads `/mnt/sdb2/lluoto/int-iScore/md_3070_0082_hi_2_3.csv`
   - reads `/mnt/sdb2/lluoto/3070_0082_repre`
   - appends results into `/mnt/sdb2/lluoto/int-iScore/md_3070_0082_hi_3_3.csv`
3. Confirmed this legacy script is **not integrated into the active package pipeline**.
4. In the active package code, frustration is still effectively a placeholder in the main modes:
   - `src/int_iscore/modes/md_snapshots.py` writes `0  # frustration_score`
   - legacy scripts also often write `frustration_score = 0`

## Current bottlenecks

### BSA

- The **active package BSA path is already on freesasa**.
- The remaining concern is not the package itself, but whether old legacy outputs were generated with
  different chain semantics (especially non-reference AF3 cases with multi-chain partner groups).

### SC

The remaining SC error is still mainly in the **surface model**, not in the final scoring formula.

Observed issues:

- The Python SC implementation can be close to CCP4 on some binary cases, but is not yet robust across
  full datasets.
- `6A6I` benchmark behavior varies strongly depending on the buried/accessibility partition strategy.
- `5HT2C` shows that treating `B+C` as a combined partner is necessary; simple pairwise averaging is not
  equivalent to CCP4 molecule-vs-molecule input.
- Multi-chain SC remains more fragile than binary SC.
- Python SES/surface generation is still an approximation, not the original CCP4 Connolly engine.

### Frustratometer integration

- `frustratometer2-master/run.py` is **not wired into the active package API**.
- The active project still lacks a clean packaged frustration-score integration path.
- Any future integration should avoid hardcoded dataset paths and should be exposed through the package mode layer rather than standalone CSV append scripts.

## Most important benchmark observations

### 5HT2C

For `pred.rank_0`:

- CCP4 `A-B`: about `0.448`
- CCP4 `A-C`: about `0.457`
- CCP4 `A-(B+C)`: about `0.450`

Representative Python values varied by implementation checkpoint, but the main conclusion remained:

- `A-B` can get reasonably close
- `A-C` is more unstable
- `A-(B+C)` must be treated as one molecule; pairwise mean is not enough

Additional geometric observation for `5-HT2C-2CD28`:

- Chains `B` and `C` are not independent distant partners; they form a connected partner-side structure.
- Direct heavy-atom minimum distance between `B` and `C` in `pred.rank_0.pdb` was approximately `0.57 Å`.
- Residue-level `B-C` contact scan found many contacts (about `99` under a `5 Å` atom-distance criterion in the inspected case).
- This strongly supports treating `B+C` as a single partner molecule for SC rather than averaging `A-B` and `A-C` independently.

### 6A6I

The full `6A6I` benchmark has been the best stress test for SC.

- Some implementation checkpoints improved representative spot cases but degraded the full dataset.
- The remaining SC problem is therefore not solved by a local heuristic alone.

## Proposed next step

### Recommended technical direction

Move the SC core to a **C++ implementation** with a true Connolly-style surface generator.

Reason:

- The current Python code still approximates the original CCP4 `sc` surface engine.
- The dominant error source is surface construction and buried/trimmed interface detection.
- A C++ implementation is the right place to reproduce:
  - convex / toroidal / concave surface points
  - CCP4-like dot assignment behavior
  - stable buried-surface partitioning
  - lower memory usage on large MD structures

### Suggested architecture

- Keep Python for orchestration and I/O
- Replace only the SC geometric core with a compiled backend
- Feed it:
  - atom coordinates
  - radii
  - chain-group selections
- Return:
  - surface-point counts
  - buried/accessibility stats
  - trimmed interface counts
  - final `Sc`

### Experimental backend note

- A pybind11 backend scaffold was integrated and builds successfully:
  - `src/int_iscore/utils/sc_backend.cpp`
- It is gated behind `INT_ISCORE_USE_SC_BACKEND=1` and is currently **disabled by default** in practice via Python fallback.
- Multiple backend iterations were tried. Build/integration is solved, but current backend geometry/scoring is still scientifically incorrect on CCP4 reference cases.
- A principled spatial-hash nearest-neighbor search was integrated into the experimental backend to replace brute-force nearest-point scoring.
- This improves backend search structure and scalability, but does **not** yet solve the main scientific mismatch; surface generation and buried/trim classification remain the dominant error sources.
- The backend therefore remains a handoff target for further work, not an active production path.

### Proposed continuation directory for another model

Use this as the active working tree:

- `/media/cuixi/data01/lluoto/int-iScore`

Most relevant next-step files:

- `src/int_iscore/utils/sc_backend.cpp`
- `src/int_iscore/utils/sc_calculator.py`
- `src/int_iscore/utils/sc_ec.py`
- `src/int_iscore/core/metrics.py`
- `calculate_md.py`
- `calculate_non_ref_af3.py`
- `calculate_af3.py`
- `frustratometer2-master/run.py`

## Files and directories involved

### Active package files

- `src/int_iscore/core/metrics.py`
- `src/int_iscore/utils/sc_calculator.py`
- `src/int_iscore/utils/sc_ec.py`
- `src/int_iscore/modes/af3_nonref.py`
- `src/int_iscore/modes/af3_reference.py`
- `src/int_iscore/modes/md_snapshots.py`

### Legacy scripts patched / inspected

- `calculate_md.py`
- `calculate_non_ref_af3.py`
- `calculate_af3.py`

### Benchmark / helper scripts

- `int_iscore_temp/run_ccp4_sc.sh`
- `int_iscore_temp/test_all_6a6i_sc.py`
- `int_iscore_temp/compare_surface_stats.py`
- `int_iscore_temp/compare_5ht2c_ccp4_python.py`
- `int_iscore_temp/analyze_6a6i_outliers.py`
- `int_iscore_temp/analyze_sc_score_terms.py`

### Frustratometer / frustration files

- `frustratometer2-master/run.py`
- `frustratometer.tar.bz2`
- `src/int_iscore/modes/md_snapshots.py` (current placeholder frustration path)

### Datasets used for comparison

- `all_3_6.csv`
- `5HT2C_intiscore.csv`
- `md_3070_0082_hi_1_3.csv`
- `md_3070_0082_hi_3_3.csv`
- `/media/cuixi/data01/lluoto/6A6I/...`
- `/media/cuixi/data01/lluoto/5-HT2C-2CD28/...`
- `/media/cuixi/data01/lluoto/split_frames/...`

## Bottom line

- **BSA**: active path is already on `freesasa`; legacy VMD path has been removed from the old scripts.
- **SC**: improved significantly from the original Python attempt, but still not a full CCP4-equivalent replacement.
- **Best next move**: reimplement the SC geometric core in C++ rather than continue stacking heuristics on the Python SES approximation.
