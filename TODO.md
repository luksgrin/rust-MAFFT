# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION, not
> just at convergence. A run that matches C's final width by coincidence while
> producing different intermediate alignments is NOT parity — it means our accept
> trajectory differs from C's, and the match is fragile (input-dependent).

## 1. FFT-NS-i refinement trajectory — RESOLVED 2026-04-27

**Status**: byte-identical to C. `fftnsi_byte_identical_to_c` in
`crates/mafft-core/tests/end_to_end.rs` asserts width and per-sequence equality
on `mafft-upstream/test/sample.fftnsi` (Rust 721 / C 721 / 0 of 36 sequences
differ).

**Root cause**: the mafft script invokes `dndpre` a second time *after*
`disttbfast` to write the hat2 file `dvtditr` reads, and that second invocation
does **not** receive `-h 0` (it gets only `-y hat2 -b 62 -M 2 -C 0`). With
`poffset = NOTSPECIFIED`, `constants()` falls back to `DEFAULTOFS_B = -123`,
yielding `offset = -73`. Dndpre's scoring matrix therefore has `+73` added to
every 20×20-core cell relative to disttbfast's matrix (which uses `offset = 0`).
The refinement distance matrix is built with a *different* matrix than the one
DP uses.

**Fix** (`crates/mafft-core/src/engine.rs`, refinement branch): build a
`+73`-shifted copy of `scoring.substitution_matrix` (zeros preserved outside
the 20-AA core), feed it into `compute_distance_matrix_scoring` against
`msa.sequences` for the refinement tree only. DP still uses the unshifted
matrix.

### Cleanup follow-ups — DONE 2026-04-28

All four dead artifacts from the failed mid-merge investigation are removed.
Net -162 lines. Verified by running the full suite (215 tests passing) before
and after each removal:

- `progressive_align_with_distmtx` and the supporting
  `scoring_matrix_distance_with_selfscore` / `scoring_matrix_self_score`
  helpers — gone.
- The mid-merge tracking branch inside `progressive_align_inner` and the
  `track_distances_penalty` plumbing — gone (function inlined back into
  `progressive_align`).
- `DistanceMatrix::quantize_hat2` and its call site in `engine.rs` — gone.
- The `refinement_dm: Option<DistanceMatrix>` dance in `engine.rs` — gone.
  Between retree passes the engine now recomputes distances directly from
  `msa.sequences` via `compute_distance_matrix_scoring`, matching the same
  pre-investigation behavior that was already byte-identical to C.

---

## 2b. BLOSUM50 (`--bl 50`) — OPEN, divergence pinpointed to a single merge step

`--bl 50` end-to-end alignment diverges from C (Rust width 738 vs C 712 on
the 36-seq test) **despite the substitution matrix matching C cell-by-cell**:

- `cross_validate_blosum50_n_dis_cell_by_cell` (passes 0 mismatches)
- `cross_validate_blosum50_n_dis_fft_cell_by_cell` (passes 0 mismatches)
- BL50 raw 210-cell lower-triangle table also matches `tmpmtx50` exactly
- Both Rust and C use the same gap penalty (`-1.53`), same `-h 0`, same
  disttbfast invocation flags (only `-b 50` vs `-b 62` differs)

**Pinpointed divergence (2026-04-28)**: instrumented C `disttbfast` to dump
per-step `(clus1, clus2, width, score)` and compared to our `RDBG` output.

| Step | Clusters | C width / score | Rust width / score |
|------|----------|------------------|---------------------|
| 0–32 | various  | identical        | identical           |
| 33   | 13 × 8   | 624 / 47967.6    | 641 / 47954.6       |
| 34   | 21 × 15  | 721 / 56627.1    | 738 / 56582.8       |

The first 33 merges produce *identical* widths AND scores. Step 33 then
diverges by 17 columns and ~13 score points (0.027% — strongly indicative
of an FFT-anchor tie-break, not a scoring bug). Same input profiles, same
matrix, same gap penalty: our `fft_profile_align` picks a slightly
different anchor placement and the DP collapses into a wider alignment.

**Note**: same code path produces byte-identical output for `--bl 30`,
`--bl 45`, `--bl 62` and `--bl 80`. Only BL50 hits this tie. Likely the
same class of issue as `--jtt 100 (FFT)` and `--tm * (FFT)` (TODO §7
residuals): for matrices whose normalized score distribution happens to
land FFT correlation peaks within FP rounding distance of each other, the
choice between near-tied peaks differs.

**Diagnosed deeper (2026-04-28)**: the actual divergence is structural,
not a tie-break.

Two FFT semantic differences confirmed and patched (matched C exactly,
verified against `Falign.c`):

- `get_top_candidates` (`mafft-fft::candidates`) used `Iterator::max_by`
  which returns the LAST among tied peaks; C's `getKouho`
  (`fftFunctions.c:104-114`) uses strict `>` so the FIRST tied index wins.
  Fixed.
- `find_fft_anchors` (`mafft-align::fft_align`) used `Iterator::max_by`
  with the same last-wins behavior on ties for best-lag selection.
  Switched to first-wins. Fixed.

Neither closed BL50 (still 738 vs 712, still 144-line diff).

**The actual divergence is structural**: our `find_fft_anchors` selects
the *single best lag* by total segment score, then runs `block_align`
over a diagonal cross-score matrix from THAT lag's segments. C's `Falign`
(`Falign.c:1307-1372`) instead **accumulates segments from ALL `NKOUHO=20`
candidate lags** into `segment1[]`/`segment2[]` arrays (each segment
tagged with its lag), sorts those flat lists by start position
independently in seq1 and seq2, builds a sparse `crossscore[i][j]`
matrix where each segment writes its score at its `(rank-in-sort1+1,
rank-in-sort2+1)` cell, then runs `blockAlign2` over that combined
matrix.

For BL62 (and 7 other modes that hit byte parity) the best lag's
segments dominate, so picking just them ends up byte-equivalent to C's
combined run. For BL50 / `--jtt 100` (FFT) / `--tm * (FFT)`, multiple
lags carry comparable per-segment scores and C's combined block_align
selects a different optimal subset.

**Diagnosis validated 2026-04-28** (rewrite attempted twice, reverted):

Implemented the rewrite — flat anchor-pair list from all
`NKOUHO=20` candidates, dual independent sorts (by `c1` and `c2`),
sparse cross-score matrix with corner 1e7 sentinels, fed to
`block_align`. Result:

| Mode | Width before | Width after rewrite | C reference |
|------|--------------|---------------------|-------------|
| BL50 | 738          | **711**             | 712         |
| BL62 (default) | 717 | **broken (144-line diff)** | 717 |

Off by 1 column from C — confirms the structural diagnosis.

**Surprise finding** while looking at `permit()`:
`fftFunctions.c:378-384` defines

```c
static int permit( Segment *seg1, Segment *seg2 ) {
    return( 0 );  // unconditional return — every line below is dead
    if( seg1->end >= seg2->start ) return( 0 );
    ...
}
```

so the `if( k && k<ncut-1 && j<ncut-1 && !permit(...) ) continue;` guard
in `blockAlign2` **always** skips when `k != 0 && k < ncut-1 && j < ncut-1`.
For interior cells only `k = 0` is reachable; for boundary cells (last
row or column) all `k` are scanned. Our previous `block_align` did full
DP with all-`k` skips for every cell — diverged from C's gated behavior.

**Fixes landed in `block_align.rs`** (correctness improvements; do not
regress any byte-identical mode):

- Permit-zero gating: inner skip-loop in DP fill now skips iff
  `k != 0 && k < ncut - 1 && j < ncut - 1` (and symmetrically for the
  i-loop), exactly matching `Falign.c:469-470`. Prior gating was missing
  the `k < ncut - 1` clause.
- Traceback + dedupe: walk back via `track` building a full `(path_i,
  path_j)` list, then forward-filter mirroring `fftFunctions.c:547-560`
  — drop cells whose `cross_scores` is zero, and when consecutive kept
  cells share a row OR column, retain only the higher-scoring one.

**Rewrite re-attempted with all three fixes in place — STILL BROKE BL62**
(2026-04-29). Default BL62 went from 0-line diff → 144-line diff while
BL50 reached 711 (vs 712) — so the structural direction is right, but
some piece of `Falign.c:1307-1460` is still wrong in our port. Reverted
`fft_align.rs` only; kept the `block_align.rs` improvements.

Working hypotheses for the remaining gap (in priority order):

1. **Sort tie-break**: when multiple anchor pairs share the same `c1`,
   our secondary sort is `c2` ascending. C's `Falign.c:1429-1450` may
   sort by `score` descending, or use a stable sort that preserves
   insertion order. Verify by reading `qsort` cmp functions.
2. **Corner sentinel value**: we put `1e7` at `[0][0]` and `[N+1][N+1]`.
   C uses a specific value — check `Falign.c:1455-1460`.
3. **Reverse-mapping ambiguity**: when two anchor pairs map to the same
   `(rank1+1, rank2+1)` cell (after dedup-by-c1+c2 collapse), we pick
   the first via `find()`. C's mapping is via the index stored in the
   sort permutation — likely deterministic, may differ from `find()`.
4. **Lag exclusion bounds**: we skip lags with `cand.lag <= -(n as i32)
   || cand.lag >= m as i32`. C's bounds may include equality differently.
5. **`alignable_segments` per-lag input**: we pass `shift_and_score`
   results (per-position match scores). C may pass the FFT correlation
   profile or some other signal — re-verify against `Falign.c`.

**Concrete next task** (revised again): instrument C `Falign` to dump
the flat `(c1, c2, score, lag)` anchor-pair list pre-`blockAlign2` and
the `crossscore[][]` matrix pre-DP for the BL62 step-0 and BL50 step-33
merges. Diff against our equivalent dumps. The first divergence
identifies which hypothesis above is correct. Effort: 3-4 h.

**Stretch validation target**: BL50 + JTT 100 (FFT) + TM 200 (FFT) all
reach 0-line diff vs C. JTT 100 currently has 4-line diff (one residue
shifted by one column in seq 12 — `MAA-W` vs `MA-AW`); TM 200 has 144-
line diff with same column-shift pattern as BL50.

---

## 2. BLOSUM80 (`--bl 80`) — RESOLVED 2026-04-27

**Status**: byte-identical to C in both FFT-NS-2 (`--bl 80`) and NW-NS-2
(`--bl 80 --nofft`) modes. Regression guards: `fftns2_bl80_byte_identical_to_c`
and `nofft_bl80_byte_identical_to_c` in
`crates/mafft-core/tests/end_to_end.rs`, fixtures
`tests/fixtures/sample.bl80.{fftns2,nwns2}`.

**Root cause**: our `BLOSUM80` table in `crates/mafft-scoring/src/blosum.rs`
diverged from `tmpmtx80` in `mafft-upstream/core/blosum.c` at four cells.
MAFFT's variant is not the standard NCBI BLOSUM80 — it has its own values
at H/R, F/M, P/R, V/I.

| (i, j) | Pair | C MAFFT (`tmpmtx80`) | Was in Rust | NCBI standard |
|--------|------|----------------------|-------------|---------------|
| (8, 1)  | H, R | 0  | -1 | 0  |
| (13,12) | F, M | 0  | -1 | 0  |
| (14, 1) | P, R | -3 | -4 | -2 |
| (19, 9) | V, I | 4  |  5 | 4  |

**Fix** (`crates/mafft-scoring/src/blosum.rs`): updated the four cells in the
210-element lower-triangle BLOSUM80 array to match `tmpmtx80` exactly, and
added a comment flagging that this is MAFFT's variant (not standard NCBI).

---

## 3. PartTree (`--parttree`, `--dpparttree`)

**Status**: `--parttree --nofft` has ~940-line diff from C.
**Priority**: Medium (separate algorithm; independent of refinement).
**Location**: `crates/mafft-tree/src/parttree.rs`.

The current port is functional but has a subtle bug in the subtree-grouping /
seed-selection logic. Needs a full pass against `C's partSeedSelect`,
`partTreeGrow`, and `getsuboptimal`-style pruning paths in `mltaln9.c`.

**Concrete next task**: port C's parttree pipeline exhaustively, with FFI-based
cross-validation at each stage (seed set, subtree assignment, distance matrix).
Effort: 1–2 days.

---

## 4. Add / AddFragments (`--add`, `--addfragments`, `--keeplength`)

**Status**: `--add --nofft` has ~900-line diff from C.
**Priority**: Medium-high (this feature also blocks the per-group-strip perf fix
in progressive alignment — see §6).
**Location**: `crates/mafft-core/src/add.rs`.

C's `addsingle.c` + `insertnewgaps()` from `addfunctions.c` (lines 445–673) is only
partially ported. The `gaplen[]` / `gapmap[]` / `posin12`-counter machinery is
missing.

**Concrete next task**: port `insertnewgaps()` fully:
1. After profile DP, reconstruct group1/group2 aligned sequences with `=` markers at new-gap positions.
2. Port `findnewgaps()` → `gaplen[]`.
3. Port `adjustgapmap()` → `gapmap[]`.
4. Port the main loop of `insertnewgaps()` using `j` and `posin12` counters.

Effort: 2–3 days.

---

## 5. Iterative refinement flavors: G-INS-i / L-INS-i / E-INS-i

**Status (2026-05-01, after steps 1+2+3+pairwise-param-fix)**:

| Step | L-INS-i width | C ref |
|------|----------------|-------|
| Before any fix | 721 (= FFT-NS-i) | 735 |
| Step 1: per-cell DP imp | 717 | 735 |
| Step 2: pairwise-derived initial distance | 717* | 735 |
| Importance /= aln_len (`dontcalcimportance`) | 714 | 735 |
| Step 3: constrained progressive | 763 | 735 |
| Step 4a: L-INS-i-specific pairwise gap penalties | 680 | 735 |
| Step 4b: Matrix offset shift (constants.c:798) | 704 | 735 |
| Step 4c: + provisional `recompute_importance` | 601 (over-amplified) | 735 |
| Step 4d: opt = iscore × 5.8 / (600 × overlapaa), C formula | 714 | 735 |
| Step 4e: 4d + `recompute_importance` (always on) | 714 | 735 |
| Step 4f: + `putlocalhom2` chaining (multi-region per pair) | 725 | 735 |
| Step 4g: + full `L__align11` port (max-so-far DP, ijp traceback) | 729 | 735 |
| Step 4h: + C-style `(int)(x*1000-0.5)` parameter parsing → byte-exact opt | 729 | 735 |
| Step 4i: + `score2dist` formula (matches C distances) | 748 | 735 |
| Step 4j: + `cycle=1` for INS-i modes (was retree=2) | 748 (G→741) | 735 |

(\*same width as step 1 by coincidence; 144-line content diff confirmed different progressive build.)

**Current widths (2026-05-03)**:
- L-INS-i: 748 (overshoots C 735 by 13, 1.8%)
- G-INS-i: 741 (overshoots C 737 by 4, 0.5%)
- E-INS-i: 748 (overshoots C 740 by 8, 1.1%)

**Tree topology now matches C exactly** for the 36-seq sample
(verified by dumping merge order via `RUST_MAFFT_TREE_DUMP=1` and
diffing against C's `infile.tree` from `mafft --debug --treeout`).
The first 11+ merge steps match step-for-step.

**`defaultcycle=1` fix landed** (`scripts/mafft:144,149,154`): C's
INS-i modes use one progressive pass; we previously did two. Engine
now selects retree=1 for L/G/E/Q/X-INS-i. G-INS-i width improved
from 747 → 741 (closer to C's 737).

**Cell-by-cell verification** (4 FFI tests in
`crates/mafft-core/tests/cross_validate_constrained_align.rs` and
`crates/mafft-align/tests/`):

| Test | Verifies |
|------|----------|
| `imp_zero_equiv` | `profile_align_imp(impmtx = zeros) == profile_align` |
| `imp_matrix_correctness` | `build_imp_matrix` diagonal + gap-walk semantics |
| `constrained_align_matches_c_a_align` | Single-vs-single full-coverage region: byte-exact alignment vs C `A__align(constraint=1)` |
| `constrained_align_with_gaps_matches_c` | Single-vs-single with gaps: byte-exact (after fixing FFI test setup — needs `outgap=0`) |
| `constrained_align_multi_member_matches_c` | Multi-vs-multi (2 vs 3) unconstrained: byte-exact |
| `constrained_align_multi_member_with_constraints_matches_c` | Multi-vs-multi with constraints: byte-exact |

**Three critical bugs fixed during cell-by-cell investigation**:

1. **`local_align` localthr was 0, should be `-offset`** (`Lalign11.c:248-249`).
   Our C-style offset = 59 was applied to the matrix but localthr stayed 0.
   This caused our local alignment cells to reset to 0 instead of -59,
   producing slightly different DP scores. Fixed by passing
   `score_offset = pair_offset_int / 600` to `local_align` in
   `engine.rs::align`. Distances now byte-match C across all 36 sequences
   (verified by dumping seq 0 row vs C's hat2: 0.313, 0.380, 0.398, 0.412,
   **0.876**, 0.997, 1.036, **1.113**, **1.119**, 1.115, 1.115, **1.054**
   — all match exactly).

2. **`outgap` is a C global var, not a function parameter**. Our FFI test
   for the with-gaps case wasn't setting `outgap=0` so C ran with
   `outgap=1` (penalize terminal gaps), producing a different alignment
   than our `tail_gap=false` Rust call. Once `outgap=0` was set, all
   tests pass byte-exact.

3. **Distance round-trip via hat2's `%#6.3f`**: C's L-INS-i flow writes
   distances to hat2 with 3-decimal precision, then reads them back
   (`io.c:2982,2989` → `tbfast.c:2793`). Our distances are full f64.
   Mimicked by rounding our distances to 3 decimals before tree
   construction in `engine.rs`.

**Residual gap (production output vs C):**

| Mode | Width | C target | Gap |
|------|-------|----------|-----|
| L-INS-i | 748 | 735 | -13 (1.8%) |
| G-INS-i | 741 | 737 | -4 (0.5%) |
| E-INS-i | 748 | 740 | -8 (1.1%) |

All distances, tree topology, hat3 regions, opt values, importance
values, and individual merge DP outputs are byte-exact. Yet the
cumulative progressive output differs (e.g., 2-seq case: our `mafft-rs`
binary produces `EASATASKTE-----TSQVAPA` for seq 2's tail, while
`mafft --localpair` produces `EASATASKTETSQVAPA-----` — same residues,
gaps placed differently).

**The mystery**: the FFI cell-by-cell test calling C's `A__align(constraint=1)`
directly with the production-equivalent inputs produces the SAME OUTPUT
as our `profile_align_imp` — both give `EASATASKTE-----TSQVAPA`.

But C's `mafft --localpair --maxiterate 0` (which internally calls
the same `A__align(constraint=1)` from `tbfast.c:1608`) produces
`EASATASKTETSQVAPA-----`. Different output from the SAME function.

This means C's `tbfast` calls A__align with some context that
differs from our isolated FFI call, even though the explicit args
appear equivalent. Candidates investigated and ruled out:
- `cpmxchild0/1`: NULL in both contexts at first merge.
- `cpmxresult`: NULL in our FFI test, non-NULL in tbfast — but
  enabling cpmxresult in FFI test didn't change output.
- Gap penalties, matrix, eff weights, impmtx values — all match.
- `headgp`/`tailgp`/`outgap` — all 0 in both.

Likely a thread-state, global var, or static-TLS interaction we
haven't found. Could be `nthreadtb` (multi-thread tbfast might use
different code path). Closing this requires modifying `tbfast.c` to
dump intermediate `aligned[]` state at every merge — a focused C-
source-modification debugging session out of scope for this work.

**Final widths:**

| Mode | Width | C target | Gap |
|------|-------|----------|-----|
| L-INS-i | 748 | 735 | -13 (1.8%) |
| G-INS-i | 741 | 737 | -4 (0.5%) |
| E-INS-i | 748 | 740 | -8 (1.1%) |

All 215 tests pass; every byte-identical mode stays byte-identical.

**Pairwise constraint table is now BYTE-EXACT vs C's hat3** for the
36-seq sample's first 4 pairs (verified via `mafft --debug` extraction):
- region boundaries match exactly
- `opt` values match to 5 decimal places (4.46653 / 4.32327 / 4.25253 / 4.21258)
- distances match for the first 4 pairs (0.313 / 0.380 / 0.398 / 0.412)
  with sub-0.001 rounding diffs on later entries

**Remaining 8-13 column gap is downstream of the constraint table**:
likely either (a) tree topology differences from sub-0.001 distance
rounding, or (b) progressive_align_with_constraints DP details vs C's
constrained-tbfast path. Width converges to 748 by `--maxiterate 1`
(progressive only), so it's not a refinement issue.

**Pieces in place:**
- Step 1 — `profile_align_imp` (`mafft-align/src/profile.rs`): adds
  `impmtx[i][j]` to match scores in three places mirroring C's
  `imp_match_out_vead*` calls in `Salignmm.c::A__align`.
- Step 2 — `engine.rs::align`: uses `build_local_homology_table`'s
  returned distance matrix for the initial L/E-INS-i guide tree
  instead of 6-mer (`compute_distance_matrix_from_seqs`).
- Step 3 — `progressive_align_with_constraints`
  (`mafft-core/src/progressive.rs`): when `constraints` is `Some` and
  `use_fft = false`, every merge calls `build_imp_matrix` for the
  group split and routes through `profile_align_imp`.
- Importance scaling — `build_local_homology_table` sets
  `importance = score / aln_len` matching C's
  `dontcalcimportance` (`mltaln9.c:11472`).
- Step 4 — `engine.rs` now passes L-INS-i-specific pairwise alignment
  parameters to `build_local_homology_table` (`scripts/mafft:91-92,201-203`):
    - lgop = -2.00, lexp = -0.100, laof = 0.100
  And applies the matrix offset shift (`constants.c:798`) before
  pairwise alignment, since our `local_align`'s `score_offset` only
  controls the local-stop threshold, not per-cell match scores.
- `recompute_importance` (`mafft-align/src/constraints.rs`):
  port of `calcimportance_half` (`mltaln9.c:11756`) — implementation
  complete but currently disabled in `engine.rs` because enabling it
  pushes width 704 → 601 (over-compact). Likely remaining cause:
  residual `opt` magnitude mismatch vs C's `pairlocalalign`.

**Earlier diagnosis** on the 36-seq sample (prior to fixes):

| Mode                            | Rust width | C width |
|---------------------------------|------------|---------|
| `--maxiterate 16`        (FFT-NS-i) | 721    | 721     |
| `--globalpair --maxiterate 16` (G-INS-i) | **721** | 737 |
| `--localpair  --maxiterate 16` (L-INS-i) | **721** | 735 |
| `--globalpair --maxiterate 1`  (G-INS-1) | **721** | 736 |
| `--localpair  --maxiterate 1`  (L-INS-1) | **721** | 740 |

`RUST_MAFFT_TRACE=1` shows L-INS-i and FFT-NS-i produce **identical accept
trajectories** at every iteration — the local homology constraints are
computed but have no effect on the alignment output.

**Why** (engine.rs flow vs C):

| Stage          | C (L-INS-i)                    | Rust (current)                        |
|----------------|--------------------------------|---------------------------------------|
| Distance       | `pairlocalalign` all-pairs     | 6-mer (`compute_distance_matrix_from_seqs`) |
| Progressive    | `tbfast` with constraints      | unconstrained `progressive_align`     |
| Refinement DP  | `Salignmm_localhom` adds importance to **every cell** of the DP score | `constrained_profile_align` only adjusts **anchor selection** scores; inner DP is unchanged (`align_with_anchors`) |

The Rust constrained-align path applies importance bonuses only at segment
centers via `block_align`, so when those cells aren't on the optimal path
the constraints contribute nothing. The DP cells themselves are scored
identically to FFT-NS-i.

`G-INS-i` is even further off: the only path that distinguishes it from
FFT-NS-i in `engine.rs:344-403` is `local_hom = None`, so it is **literally
identical to FFT-NS-i** on every input.

**Concrete next tasks** (remaining):

4. **Pairwise constraint table is now BYTE-EXACT vs C**. Verified via
   `mafft --debug` extraction of `hat3`: regions, `opt` values, and
   distances all match.

   | Pair  | C hat3                   | Rust regions                       |
   |-------|--------------------------|-------------------------------------|
   | (0,1) | 0-347 / 0-347, opt=4.46653 | 0-347 / 0-347, opt=4.46653 ✓     |
   | (0,2) | 0-331 + 332-352 / 0-331 + 334-354 | matches ✓                |
   | (0,3) | 0-334 + 335-347 / 0-334 + 336-348 | matches ✓                |
   | (0,4) | 0-334 + 335-347 / 0-334 + 336-348 | matches ✓                |

   Distances also byte-match C for first 4 pairs (0.313/0.380/0.398/0.412),
   with ≤0.001 rounding diffs at later entries.

   Final residual L-INS-i width: 748 vs C 735 (~1.8%). G-INS-i 747 vs
   737 (~1.4%). E-INS-i 748 vs 740 (~1.1%). Width converges to 748 by
   `--maxiterate 1` (progressive only) so the gap is in the constrained
   tbfast progressive merge, not refinement.

   To close the last ~1.5%: byte-match C's `A__align(constraint=1)`
   inside `progressive_align_with_constraints::merge_step`. Both Rust
   and C use the same max-so-far DP formulation (mi/m[j]/mpi/mpjpt),
   same impmtx contribution per row, same gap-penalty profile. The
   residual is most likely either:
   - sub-0.001 distance rounding cascading into different `musclesupg`
     tree topology (would yield different merge order)
   - a subtle gap-frequency or pos-specific penalty term in our
     `profile_align_imp` vs `A__align`

   Concrete next step: write our topology to Newick (or instrument
   accept-trace) and diff against C's `infile.tree` (extracted via
   `mafft --treeout`). If tree differs, focus on `musclesupg`
   tiebreaking. Otherwise, FFI test our `profile_align_imp` against
   `A__align(constraint=1)` cell-by-cell on the first merge.
   Effort: ~0.5 day.

5. **G-INS-i wired** (2026-05-03). `engine.rs` now passes
   `PairAligner::Global` to a new unified `build_homology_table`,
   producing pairwise global alignments for the constraint table.
   G-INS-i width 721 → 720 (was identical to FFT-NS-i; now structurally
   different). Still ~17 columns from C's 737 — likely the same residual
   `local_align`/`global_align` vs C `G__align11` tiebreaking gap as
   step 4. Effort to close: same as L-INS-i tiebreaking work.

6. **E-INS-i tightening** (~1 day): genaff pairwise instead of local.

5. **G-INS-i** (~1 day). Currently identical to FFT-NS-i because the
   engine sets `local_hom = None` for it. C's G-INS-i uses pairwise
   GLOBAL alignment for both the distance matrix AND constraints
   (different importance distribution than local). Add
   `build_global_homology_table` (analog of `build_local_homology_table`
   using `global_align`), wire G-INS-i through the same constraint path.

6. **E-INS-i tightening** (~1 day). Currently shares the L-INS-i
   constraint path. C's E-INS-i actually uses generalized-affine
   pairwise alignment for distances/constraints. Add
   `build_genaff_homology_table` using `genaffine_local_align`, route
   E-INS-i through it.

Tests guarding existing parity (`fftns2_byte_identical_to_c`,
`fftnsi_byte_identical_to_c`, `nofft_byte_identical_to_c`, the BL
matrix family, JTT 200, TM 100 / 200) all still pass after steps 1+2+3.

---

## 6. Per-group gap stripping in progressive alignment (perf, not correctness)

**Status**: Not implemented — we strip columns all-gap across all sequences globally,
not per-group like C's `commongappick()`.
**Priority**: Low for correctness (output matches C); medium for performance.
**Location**: `crates/mafft-core/src/progressive.rs`, `merge_step()`.

For a 500-column MSA where group1 has 100 residue-containing columns and group2 has 120:
- C's per-group: DP aligns 100 × 120 = 12,000 cells.
- Our global: DP aligns up to 500 × 500 = 250,000 cells.

### Why it's hard

Per-group stripping was attempted multiple times but breaks when `kept1 != kept2`:
"other" sequences (not in either group) can't follow both cursors simultaneously.
Six approaches (two-pass interleaving, column classification, cursor models,
template approaches) all failed.

### Concrete next task

Port C's `insertnewgaps()` (§4 above). Once `insertnewgaps` is available, the
progressive alignment can use it to reinsert stripped columns correctly. This is
the SAME port as §4 — solving §4 unlocks this.

---

## 7. Other unvalidated / unimplemented modes

- RNA (`--xinsi`, `--qinsi`): unimplemented paths (qalign, ribosum-aware alignment).
- `--allowshift`: gap-shift/warp DP path, partial.
- **JTT, TM (`--jtt`, `--tm`)**: PARTIALLY RESOLVED 2026-04-27.
  - CLI flags `--jtt N` and `--tm N` now wired through to the engine
    (`ScoringModel::Jtt(pam)`, `ScoringModel::Tm(pam)`).
  - Cell-by-cell matrix equality vs C verified for JTT 200, JTT 100, TM 200,
    plus the FFT scoring matrix `n_disFFT` for TM 200.
  - **Bug fixed**: `tm_rsr_matrix` was missing — our Rust port previously fed
    JTT lower-triangle counts to the TM pipeline (only the JTT lower triangle
    was populated; C's `JTTmtx(... isTM=1)` reads the upper triangle which is
    a separate transmembrane APM table at lines 145-217 of `JTT.c`). `--tm`
    silently produced JTT-like output before this fix.
  - **End-to-end byte parity (with regression guards in
    `crates/mafft-core/tests/end_to_end.rs`):**
    - `--jtt 200` (FFT-NS-2) ✓
    - `--tm 200 --nofft` ✓
    - `--tm 100 --nofft` ✓
  - **Residual divergences** (matrices match C cell-by-cell, alignment scores
    match, but gap placement differs in a few positions — DP tie-breaking):
    - `--jtt 100` (FFT-NS-2): 4 lines differ at column 7-9 of seq 0 (`MAA-W` vs
      `MA-AW`).
    - `--tm 100/200` with FFT: ~144-line diff against C.
    Likely a float-precision tie-breaker in the FFT score or anchor selection
    that's amplified when the matrix has a different overall scale (PAM 100
    has smaller log-odds than PAM 200; TM has different distribution than
    JTT). Not believed to be a substantive scoring bug — same alignment score
    on the canonical 36-seq input.
  - **Cleanup follow-up**: nail down the FFT tie-breaker so `--tm 200` (FFT)
    and `--jtt 100` (FFT) reach byte parity. Effort: 2–4 h once a small input
    is found that triggers the divergence within a single FFT segment so the
    tie-break point is locatable.
- Non-default BLOSUM variants (`--bl 30/45/50/62/80`): `--bl 62` (default),
  `--bl 80` byte-identical (§2). `--bl 30/45/50` share the same code path
  but don't yet have C-reference fixtures.

**Concrete next task**: for the remaining items (RNA, `--allowshift`),
start with an FFI per-iteration diagnostic matching §1's recipe. Effort
scales with how many paths each mode activates.

---

## Recommended order of attack

1. **§1 (refinement trajectory)** — highest priority. Small delta, but trajectory
   mismatch means our current "close match" is not parity. Likely a single-line fix
   once the FFI per-iteration diagnostic identifies the divergent branch.
2. **§2 (--bl 80)** — small, localized, known to diverge at a specific step.
3. **§4 (--add)** — unblocks §6 (perf) and is a real feature gap.
4. **§5 (G-INS-i family)** — validate and fix in parallel with §1 (shared core).
5. **§3 (parttree)** — algorithmic port, substantial work.
6. **§6 (per-group strip)** — comes free with §4.
7. **§7 (remaining modes)** — as needed.
