# FastSAS Optimization Challenge: Fix Negative Correlation with CCP4 SC

## Context

I have a C++17 high-throughput Shape Complementarity (SC) engine with two modes:

1. **FastSAS** (mode 0): Solvent-Accessible Surface + KD-tree scoring. ~250 ms/complex.  
   Designed as a coarse filter for ranking 10,000+ protein complexes.

2. **Accurate** (mode 1): Calls Python SCASA bridge via popen, which computes CCP4-compliant Connolly molecular surface SC. ~10 s/complex.  
   This is the gold standard, matching CCP4 within +/-2%.

**Goal**: FastSAS should rank complexes with at least *positive* Pearson correlation with Accurate (CCP4) SC scores.

## The Problem

FastSAS is **negatively correlated** with Accurate SC across 4 AlphaFold3-predicted protein complexes (6gs2, 6h4b, 6ii6, 6ono), each with 1000 seed variations. After multiple parameter tuning attempts, the correlation has only gotten worse:

| Trial | Weight (w) | Cutoff | Interface Dist | Aggregation | Pearson r |
|-------|-----------|--------|----------------|-------------|-----------|
| Original | 0.08 | 3.0 A | 8.0 A | median | -0.416 |
| v2 | 0.08 | 3.0 A | 8.0 A | mean | -0.452 |
| v3 | 0.20 | 3.0 A | 5.0 A | mean | -0.476 |
| Gaussian | 0.05 | 4.0 A | 8.0 A | median | -0.382 |

The Gaussian Overlap trial (commit 23a9b8c) changed the formula to
`normal_score * exp(-w * d^2)` where `normal_score = max(0, -n_dot)`, making
normal-vector alignment the dominant factor. This is the best effort to date,
improving correlation by 8% from baseline, but remains negative.

This means when FastSAS says a complex has high SC, the CCP4-correct score tends to be low, and vice versa. The ranking is completely inverted.

## Detailed Data

### Per-PDB Statistics (1000 samples each for FastSAS original)

| PDB ID | FastSAS Mean | FastSAS Std | Accurate Mean | Rank(Fast) | Rank(Acc) |
|--------|-------------|------------|---------------|------------|-----------|
| 6gs2 | 0.409 | 0.058 | 0.570 | #1 | #4 |
| 6h4b | 0.341 | 0.053 | 0.577 | #3 | #3 |
| 6ii6 | 0.377 | 0.056 | 0.587 | #2 | #2 |
| 6ono | 0.323 | 0.053 | 0.650 | #4 | #1 |

The per-PDB ranking is almost perfectly inverted. Overall paired data (N=80): Pearson r=-0.416, FastSAS mean=0.361, Accurate mean=0.596.

### Paired Sample Detail (10 of 80)

```
6gs2_seed-1020037_sample-0: Fast=0.459  Acc=0.585  ratio=78%
6gs2_seed-1020037_sample-1: Fast=0.307  Acc=0.570  ratio=54%
6gs2_seed-1020037_sample-2: Fast=0.246  Acc=0.595  ratio=41%
6gs2_seed-1020037_sample-3: Fast=0.407  Acc=0.589  ratio=69%
6gs2_seed-1020037_sample-4: Fast=0.285  Acc=0.593  ratio=48%
6gs2_seed-1042502_sample-0: Fast=0.482  Acc=0.584  ratio=82%
6gs2_seed-1042502_sample-1: Fast=0.440  Acc=0.560  ratio=79%
6gs2_seed-1042502_sample-2: Fast=0.351  Acc=0.576  ratio=61%
6gs2_seed-1042502_sample-3: Fast=0.457  Acc=0.577  ratio=79%
6gs2_seed-1042502_sample-4: Fast=0.318  Acc=0.562  ratio=57%
```

Key observation: within the same PDB (e.g., all 6gs2 seeds), Accurate SC is tightly clustered (0.56-0.60), while FastSAS varies widely (0.25-0.48). This suggests FastSAS is capturing seed-level structural noise rather than true interface quality.

## Code: FastSAS Scoring Function

The core scoring function (pseudocode):

```
// 1. Generate SAS surface dots for each chain pair
//    Fibonacci sphere sampling on each atom surface + 1.4A probe radius
//    Occlusion: remove dots within neighbor atom VdW radius
fast_sas_surface(iface_atoms, self_kdtree, all_atoms, dots);

// 2. Classify dots: buried vs accessible
classify_buried_fast(dots, partner_kdtree, buried, accessible);

// 3. Score: for each buried dot in chain1, find NN in chain2 buried dots
//    score = aggregate( -dot(n1,n2) * exp(-w * d_squared) )
//    SC = (score_ab + score_ba) / 2

for each dot i in buried1:
    nn_idx, dist_sq = kdtree2.knnSearch(dot_i, 1)
    if sqrt(dist_sq) > CUTOFF: skip
    s_val = -dot(n1_i, n2_nn) * exp(-WEIGHT * dist_sq)
    scores.push(s_val)
return aggregate(scores)
```

**Key parameters:**
- PROBE_RADIUS = 1.4 A
- INTERFACE_DISTANCE = 8.0 A (atoms within this of partner are "interface")
- DOT_DENSITY = 15 dots/A^2 (Fibonacci sphere sampling)
- WEIGHT_FAST = 0.08 (weight in exp(-w * d^2))
- FAST_DISTANCE_CUTOFF = 3.0 A (hard truncation)

**Surface difference:** FastSAS uses Solvent-Accessible Surface (atom center + VdW radius + probe radius). Accurate uses Connolly Molecular Surface (probe rolling over VdW spheres). SAS dots are approximately 1.4 A farther from atoms than Connolly dots, leading to inflated NN distances.

## Accurate (Reference) Scoring

The Accurate mode calls Python SCASA, which computes the full Connolly mds (convex + toroidal + concave patches) with analytical normals, then applies the CCP4 formula:

```
S = sum( -dot(n1, n2) * exp(-0.5 * d^2) ) / N_ref_dots
```

Parameters: w=0.5, trim=1.5 A, interface=2.0 A. This matches CCP4 within +/-2%.

## What We Have Tried

1. **Weight (w)**: 0.08 to 0.5 — results near zero (SAS distances too large)
2. **Weight (w)**: 0.08 to 0.2 — correlation worsened
3. **Aggregation**: median to mean — correlation slightly worsened
4. **Interface distance**: 8.0 to 5.0 — fewer dots, correlation worsened
5. **Cutoff removal**: 3.0 to none — all scores collapsed near zero

5. **Gaussian Overlap** (Phase 4, commit 23a9b8c): Changed formula to
   `normal_score * exp(-w * d^2)` with `normal_score = max(0, -n_dot)`, w=0.05,
   cutoff=4.0 A. Result: r=-0.382 (+8% from baseline). This is the best attempt --
   normal-vector-dominant scoring partly compensates for SAS distance inversion
   but cannot fully overcome it.

## Physical Diagnosis

The root cause is **SAS surface distance inversion**:

- SAS dots are ~1.4 A farther from atoms than Connolly dots
  (atom center + VdW radius + 1.4 A probe radius)
- At perfect complementarity, SAS surfaces are pushed apart, causing NN distances
  to increase, which gives **lower** scores in the exp(-w * d^2) weighting
- This is the opposite of Connolly SC, which gives higher scores at tighter packing
- No parameter tuning within the SAS surface framework can fix this -- it requires
  a fundamentally different surface type (e.g., convex-only Connolly approximation)

## What We Need Help With

**Primary question**: Can FastSAS (SAS-based surface) be made even *positively* correlated with CCP4 SC (Connolly-based surface) through algorithmic changes?

**Specific areas to analyze:**

1. **Surface type**: Is SAS fundamentally incapable of capturing SC? Would a faster Connolly approximation (e.g., convex-only dots + radial normals) work better?

2. **Scoring formula**: Is -dot(n1,n2)*exp(-w*d^2) the wrong formula for SAS surfaces? Should we use:
   - A different distance metric (atom-to-atom instead of dot-to-dot)?
   - A different weight function (step function, linear decay)?
   - Normalization by interface area or dot count?

3. **Buried dot selection**: The current classify_buried_fast checks if a SAS dot is within the partner atom VdW radius. This might misclassify interface vs non-interface dots.

4. **Within-PDB variance**: FastSAS varies 2x within the same PDB while Accurate is stable. How can we make it more robust?

6. **Feature-based approach**: Instead of direct scoring, could we extract features from the SAS surface (interface area, dot density, curvature) and train a simple regressor against Accurate SC?

## Constraints

- Must remain C++17, single-file (sc_hts.cpp)
- Must remain under 300 ms/complex (cannot do full Connolly mds)
- Must use nanoflann KD-tree (already integrated)
- Output CSV to stdout

## Relevant Files (on remote server)

- sc_hts.cpp — main C++ source (~1000 lines)
- src/int_iscore/utils/sc_calculator.py — Python SCASA reference
- src/int_iscore/utils/connolly.py — Connolly mds implementation
- sc_bridge.py — Python bridge for Accurate mode
- /tmp/fast_results.csv — 4,614 FastSAS results (original params)
- /tmp/acc_results.csv — 80 Accurate results for comparison
- /tmp/fast_v3.csv — 80 FastSAS v3 results
- pdb_output/ — 12,000 PDB files (12 PDBs x 1000 seeds), all A-B pairs

## Test Command (on remote server: cuixi@10.19.25.48)

```bash
cd /media/cuixi/data01/lluoto/int-iScore
g++ -std=c++17 -O3 -fopenmp -I. sc_hts.cpp -o sc_hts
./sc_hts /media/cuixi/data01/lluoto/6A6I/seed-36904_sample-0/model.pdb A B 0 test
```
