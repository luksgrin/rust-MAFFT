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

### Cleanup follow-ups (low priority, none affect parity)

- `progressive_align_with_distmtx` and the supporting
  `scoring_matrix_distance_with_selfscore` /
  `scoring_matrix_self_score` helpers were added to track mid-merge distances
  before the real cause was found. They're still wired through `engine.rs`'s
  retree loop but their result is now discarded (`refinement_dm.take()` is a
  no-op). The whole mid-merge tracking path can be removed.
- The `quantize_hat2` call on the refinement DM is harmless but unnecessary
  with the offset-shift fix in place — dndpre writes hat2 with `%5.3f`, but
  dvtditr's tree builder (`fixed_musclesupg_double_treeout`) reads them back
  with full atof precision; the 3-decimal round-trip is encoded in the file
  format, not in the in-memory pipeline we're emulating. Worth verifying with
  a one-line ablation that removing `dm.quantize_hat2()` keeps the test green.

---

## 2. BLOSUM80 (`--bl 80`)

**Status**: diverges from C under `--nofft` (width 700 vs C's 712). Raw BLOSUM80
data and normalization are correct (shared code path with BLOSUM62, which matches
C exactly).
**Priority**: Medium.
**Location**: `crates/mafft-scoring/src/blosum.rs`, `crates/mafft-core/src/progressive.rs`.

Divergence starts at retree 1 step 19. Likely cause: an ambiguity-code (B/Z/X)
cell that differs between BLOSUM80 and BLOSUM62 expansions, hitting at a sequence
position where BLOSUM62 sequences happen to have a residue.

**Concrete next task**: add a per-retree-step RDBG comparison for `--bl 80 --nofft`
that diffs profile/gap-counts/alignment output at each step. The first divergent
step will pinpoint the matrix or scoring cell responsible. Effort: 2–3 hours.

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
- JTT, TM (`--jtt`, `--tm`): matrix construction FFI-validated cell-by-cell; alignment trajectory unvalidated.
- Non-default BLOSUM variants (`--bl 30/45/50/62/80`): share codepath with BLOSUM62; `--bl 80` diverges (§2); others unverified.

**Concrete next task**: for each, start with an FFI per-iteration diagnostic
matching §1's recipe. Effort scales with how many paths each mode activates.

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
