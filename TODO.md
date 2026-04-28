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

## 2b. BLOSUM50 (`--bl 50`) — OPEN MYSTERY 2026-04-28

`--bl 50` end-to-end alignment diverges from C (Rust width 738 vs C 712 on
the 36-seq test) **despite the substitution matrix matching C cell-by-cell**:

- `cross_validate_blosum50_n_dis_cell_by_cell` (passes 0 mismatches)
- `cross_validate_blosum50_n_dis_fft_cell_by_cell` (passes 0 mismatches)
- BL50 raw 210-cell lower-triangle table also matches `tmpmtx50` exactly
- Both Rust and C use the same gap penalty (`-1.53`), same `-h 0`, same
  disttbfast invocation flags (only `-b 50` vs `-b 62` differs)

Yet the alignment DP picks different equivalent paths. The scoring inputs
are provably identical, so the divergence has to be either (a) a different
hidden global that varies with `nblosum` (unidentified), or (b) the C DP
makes a different tie-breaking decision under BL50's specific number
distribution.

**Note**: same code path produces byte-identical output for `--bl 30`,
`--bl 45`, `--bl 62` and `--bl 80`. Only BL50 diverges.

**Concrete next task**: instrument C's `MSalignmm` (or whichever DP path
handles BL50 here) with per-cell logging at the first divergent merge step,
compare to Rust's per-cell DP scoring. Effort: 4-6 h.

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

**Status**: implemented but not cross-validated against C per-iteration.
**Priority**: Medium. Likely in a similar state to FFT-NS-i (close but not iteration-exact).
**Location**: `crates/mafft-core/src/engine.rs` (mode selection), `crates/mafft-core/src/refinement.rs` (shared pipeline).

**Concrete next task**: run one canonical example per mode against C and diff
per-iteration widths, as in §1. The fixes may be shared with §1 (same refinement
core) or may surface mode-specific paths (e.g., constrained DP for L-INS-i).
Effort: 0.5–1 day per mode.

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
