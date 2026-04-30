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

**Status (2026-04-30)**: ALL THREE MODES PRODUCE FFT-NS-i OUTPUT. Diagnosed
on the 36-seq sample:

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

**Concrete next task** — three independent ports, in size order:

1. **Constraint-aware refinement DP** (~1 day): port C's
   `Salignmm_localhom` so per-cell match scores in the inner DP add the
   homology importance. Test target: L-INS-i diverges from FFT-NS-i.
2. **Distance from pairwise alignments** (~0.5 day): port the
   `pairlocalalign` distance loop (or reuse `local_align` results from
   `build_local_homology_table` to derive per-pair distance =
   `1 - identity`). Test target: L-INS-i guide tree differs from
   FFT-NS-i's. Already partially in `build_local_homology_table` —
   rewire `engine.rs` to use it as the initial DM for INS-i modes.
3. **G-INS-i / E-INS-i variants** (~1 day each): pairwise GLOBAL
   alignment for G-INS-i (use existing `profile_align`/`global_align`),
   pairwise GENAFF for E-INS-i (needs a generalized-affine DP variant
   not currently in `mafft-align`).

Steps 1-2 close L-INS-i; step 3 adds G/E variants. Total effort ~3-4
days for byte-parity.

Mid-step validation: even before byte parity, after step 1 the Rust
L-INS-i width should *change* relative to FFT-NS-i — that's the first
visible signal that constraints are flowing through the DP.

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
