# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION, not
> just at convergence. A run that matches C's final width by coincidence while
> producing different intermediate alignments is NOT parity — it means our accept
> trajectory differs from C's, and the match is fragile (input-dependent).

## 1. Refinement alignment trajectory (Rust diverges from C at every iteration)

**Status**: Rust final=720, C final=721 (0.1% gap) — but this is coincidental. The
per-iteration widths diverge from C at every iteration:

| maxiterate | Rust | C   | delta |
|------------|------|-----|-------|
| 1          | 717  | 721 | -4    |
| 2          | 716  | 719 | -3    |
| 3          | 723  | 722 | +1    |
| 4          | 720  | 721 | -1    |
| 5+         | 720  | 721 | -1    |

**Priority**: HIGH — the trajectory mismatch indicates our accept decisions at some
branches diverge from C's, meaning subsequent branches run on different state than C
does. We happen to converge near C's final width, but on a different input this
could diverge more.

**Location**: `crates/mafft-core/src/refinement.rs` — `realign_all()` and `iterative_refine()`.

### What is already confirmed matching C (FFI cross-validation)

Every individually-testable component is byte-identical to C:
1. `BranchWeights` / `weightFromABranch` — exact to 2.22e-16 across all 70 branches.
2. Per-group normalization / `fastconjuction_noname` — exact.
3. Profile construction / `cpmx_calc_new` + `st_*GapCount` + `gapcountf` — exact.
4. Intergroup scoring / `intergroup_score` — exact.
5. `alignable_segments` = C's `alignableReagion(lag=0)` — 0/69 divergent centers across all iter-0 branches.
6. Full segmented Falign pipeline (`alignable_segments` + per-segment `commongappick` + `profile_align`) — 0/69 divergent outputs across all iter-0 branches, both clean and evolving state (FFI tests `compare_segmented_falign_with_c_iter0` and `evolving_state_rust_vs_c_falign_iter0`).

### Fixes landed this session

- **MINLEN clamp** in `BranchWeights::new` matching `checkMinusLength` (0.001).
- **Nongap_freq boundary** in `profile_align`: `nongap_freq[n] = 1.0` matching C's `gapfreq1pt[lgth1]` padding.
- **`compute_split_score` sequential sum** (no `par_iter`) for deterministic accumulation.
- **`realign_all` rewritten to match `dvtditr` kobetsubunkatsu=1 flow**: no FFT, `alignableReagion(lag=0)` for cut points, per-segment `commongappick` + `profile_align`.
- **Boundary gap state (`sgap`/`egap`)** corrections on `ogcp[0]` and `fgcp[length-1]` per segment matching `new_OpeningGapCount`/`new_FinalGapCount`.
- **Identity check uses 2 representatives only** (`s1 = group1[0]`, `s2 = group2[0]`) matching `tditeration.c:2184-2185`. Our previous all-sequences check produced spurious accepts when column rearrangements preserved the representatives — this single fix closed 9 of 10 columns (729 → 720).

### What the per-iteration FFI diagnostic discovered (this session)

A comprehensive diagnostic (`trajectory_rust_vs_c_all_iters` and
`real_rust_refine_vs_c_falign_iter0_bytediff` in
`crates/mafft-core/tests/trace_refinement.rs`) ran Rust's real `iterative_refine`
branch-by-branch alongside C's `Falign` via FFI with identical weights, inputs
and scoring. Findings:

- **For every single branch of iter 0, Rust's output matches C's `Falign` output
  byte-identically** (`bytes_match_c=true` at all 9 accepts). Identity decisions,
  `old_score`, and `tscore` computed through `intergroup_score` FFI all match
  to machine epsilon across both implementations.
- **The sequence of accept/reject decisions matches C exactly**: same 9 accepts
  at the same branches with the same score deltas.
- **`fixed_musclesupg_double_treeout` topology matches our `musclesupg`
  byte-for-byte** (0/35 mismatches, verified in `tree_topology_rust_vs_c`): the
  tree C builds from the hat2 distance matrix is identical in both cluster
  membership and branch lengths to the tree Rust builds from the same matrix.
- **Yet the final state after iter 0 still differs from real `dvtditr -I 1` by
  4 columns**: Rust converges to 717, dvtditr-via-mafft-script converges to 721.

Direct `dvtditr` invocation (bypassing the mafft script) was used to bisect:

| Command | Width | Notes |
|---------|-------|-------|
| `dvtditr ... -b 62 -I 1 ...` (full mafft args, -C 1 BAATARI2) | 721 | C real behavior |
| `dvtditr ... -b 62 -I 1 -p BAATARI1 ...` | 721 | serial path, same result |
| `dvtditr ... -b 62 -I 1 -C 0 ...` (non-threaded path) | 721 | same result |
| `dvtditr ... -b 62 -I 1 -C 8 ...` (8-thread BAATARI2) | 728 | threads differ |
| Our FFI `Falign` per-branch simulation with niter=1 | 717 | matches Rust real refine |

So the difference is NOT in:
- `alignableReagion` (segment centers match exactly)
- `Falign` per-segment (outputs match byte-for-byte via FFI)
- `MSalignmm` / `A__align` / profile DP
- Weights (`weightFromABranch` matches to 2.22e-16)
- Distance matrix (`dndpre` formula matches byte-for-byte)
- Tree topology (`fixed_musclesupg_double_treeout` matches 0/35)
- `intergroup_score` (matches to machine epsilon)
- The identity check (2-reps)
- The accept threshold logic

Yet the `dvtditr` binary produces different output than the sum of its parts
called via FFI. The residual must be in something the `dvtditr` binary does
AROUND or BEFORE the per-branch loop that our FFI simulation is missing —
likely a subtle state-management detail in `TreeDependentIteration` itself.

### Remaining suspects (narrowed by the diagnostic)

1. **`commongappick(njob, seq_g)` at dvtditr.c:717 on the full alignment BEFORE
   iteration starts.** Verified no-op for our sample (no all-gap columns). Matters
   for other inputs but not this one.
2. **`searchAnchors` at dvtditr.c:891 top-level segmentation.** Returns 0 anchors
   for our sample (`c_search_anchors_on_sample` test) → whole alignment is one
   outer segment → matches our behavior. Matters for other inputs.
3. **Threading path.** Non-threaded `-C 0` still gives 721; `-C 1 -p BAATARI2`
   uses the multi-thread path (with 1 worker) but also gives 721. `-C 8`
   gives 728. So non-determinism exists, but all single-worker paths (threaded
   or not) converge to 721 — which is what we fail to reach.
4. **`TreeDependentIteration`'s per-branch work beyond Falign.** This function
   does more than call Falign — it manages `aseq`/`bseq` state, calls `fastconjuction_noname`
   (our `rw1n`/`rw2n` building), does score tracking with `history[iter][l][k]`,
   and updates `converged`. Something in this management differs from our
   Rust `iterative_refine`, even though every individual operation matches.
5. **Multi-threaded path internals.** When dvtditr uses the `#ifdef enablemultithread`
   path with 1 worker, the `athread` function has its own branch-queueing, gain
   tracking, and acceptance protocol (BAATARI2 strategy). This is a DIFFERENT
   algorithm from the sequential branch-by-branch loop, not just a parallel
   implementation of it.

### Investigations that did NOT fix the divergence

- **hat2 3-decimal precision round-trip** (`DistanceMatrix::quantize_hat2`). C's
  `WriteFloatHat2_pointer_halfmtx` stores distances with `DFORMAT = "%#6.3f"` and
  `dvtditr` reads them back via `atof`, quantizing the refinement distance matrix
  to 3 decimals. We now apply the same quantization in `engine.rs` before building
  the refinement tree (matches C's behavior exactly). This did NOT change the
  final width (still 720 vs C's 721) — quantization was a consistency fix, not the
  trajectory fix.
- **Weight computation pipeline.** Rebuilding weights via our own topology and via
  C's `weightFromABranch` FFI with our topology gives identical results to
  machine epsilon; the tree topology itself matches C byte-for-byte *when fed the
  same distance matrix on both sides*.
- **athread BAATARI2 branch ordering.** With `-C 0 -p BAATARI2` (the flag
  combination the mafft script actually passes for `--maxiterate 100`) C takes the
  sequential non-threaded path in `TreeDependentIteration`, which alternates
  l-direction per iteration (tditeration.c:1641-1648). This matches our Rust code
  exactly. The multi-thread path (`nthread > 0`) does differ — `branchtable` walks
  always forward — but that path isn't exercised under `-C 0`.

### Root cause (confirmed this session)

C's refinement tree is built from distances computed DURING the progressive
merge, not from the final progressive alignment. `disttbfast.c` line 2151:
```
newdistmtx[i][j] = (1.0 - naivepairscorefast(mseq1[i], mseq2[m], ..., penalty_dist) / bunbo) * 2.0
```
where `mseq1[i]` and `mseq2[m]` are the CURRENT progressively-aligned profile
sequences at the step that first merges i and m. Every pair gets its distance
recorded exactly once — at the moment the two sequences' clusters first meet.

Our Rust computes distances from the FINAL alignment (after all merges), so we
emit systematically different numbers. Diagnostic `refinement_distance_vs_c_hat2`
in `tests/end_to_end.rs` (fixture: `sample.fftnsi.hat2` extracted from C's
`--debug` run) shows:
- Rust d(0,1) quantized = 0.277, C hat2 d(0,1) = 0.248
- 595/630 pairs differ by ≥0.05
- max |rust_q - c_hat2| = 0.174

The refinement tree therefore differs in branch lengths, which perturbs per-branch
weights (`weightFromABranch` is length-driven), which perturbs Falign's DP choices.
Our tree is a *better* tree in a sense (uses full-precision final distances), but
we're chasing byte parity, not optimality.

### Status after mid-merge distance port (2026-04-24)

Mid-merge distance tracking **is implemented** (`progressive_align_with_distmtx`
in `crates/mafft-core/src/progressive.rs`, plumbed through `engine.rs` so every
retree pass emits a `DistanceMatrix` for the next pass's tree, and the final
pass's matrix is fed to `musclesupg` for the refinement tree after
`quantize_hat2`). Supporting helpers
`scoring_matrix_distance_with_selfscore` and `scoring_matrix_self_score` live in
`mafft-tree::distance`. Diagnostic `midmerge_distmtx_matches_c_hat2` in
`tests/end_to_end.rs` compares our mid-merge matrix to C's `sample.fftnsi.hat2`.

**Result**: trajectory is unchanged — Rust fftnsi still 720, C 721. Our
mid-merge matrix still diverges from C's hat2: 627/630 pairs with |diff| ≥
0.005, max |diff| = 0.174, our d(0,1) = 0.277 vs C 0.248, d(29,30) = 0.011 vs
C 0.010. So: infrastructure is correct, the formula matches C's 2151 line
exactly — but the actual *sequences* we feed into that formula at the joining
step differ from C's, even though our final progressive output is byte-identical
to C's (`fftns2_byte_identical_to_c` still passes).

Concretely: at the step that first joins sequences 0 and 1 (our step 13, both
groups at width 353), the profile rows we have for seq 0 and seq 1 score
differently from C's rows at the corresponding step, despite both implementations
converging to the same final alignment. Either the merge order differs (same
final alignment via a different tree), or the per-step DP choices differ in
ways that cancel by the end.

### New concrete next task

**Find where our mid-merge state diverges from C's.** Candidates:
1. Compare our pass-1 (initial) topology against C's via FFI — if trees differ
   but both produce the same final alignment, that's the signal. Instrument
   `dump_rust_topology`-style printouts side-by-side.
2. If topologies match, instrument `merge_step_cached` to dump aligned\[group1[0]\]
   after each step, and a matching `fprintf` in disttbfast.c's inner loop after
   each Falign, then diff the first step where they diverge.
3. A single pass of `RUST_MAFFT_TRACE_DM=1` (handle in `progressive.rs`) already
   prints the mid-merge distances by step — extend to also print merge group
   contents for rapid cross-check.

### Previous concrete next task (now done, did NOT close the gap)

**Port disttbfast's mid-merge distance tracking to our progressive aligner.** At
every merge step, after DP produces the joined profile, compute
`d(i,j) = (1 − naivepairscore11(mseq1[i], mseq2[j], penalty_dist) / min(self_i, self_j)) * 2.0`
for each (i in group1) × (j in group2), and record it in a distance matrix indexed
by the original sequence IDs. That matrix (quantized via `quantize_hat2`) becomes
the refinement tree input, replacing the current `compute_distance_matrix_scoring`
call in `engine.rs`. Hooks needed:
1. `progressive.rs::merge_step` returns (or populates) `(i, j, d)` records for the
   pair it just merged.
2. `engine.rs` collects those records across all merges into a `DistanceMatrix`.
3. That matrix, not `compute_distance_matrix_scoring(&msa.sequences, ...)`, feeds
   `musclesupg` for the refinement tree.
4. Keep `quantize_hat2()` applied after.

Effort: 3–6 hours. Validation: `refinement_distance_vs_c_hat2` drops to ~0 diff,
`fftnsi_width_matches_c` matches at 721 and all 36 sequences identical.

### Alternative concrete next task (smaller scope, different branch)

**Read tditeration.c multi-thread path (athread) end-to-end** and understand
what BAATARI2 actually does. The single-worker multi-thread path (`-C 1 -p BAATARI2`)
produces 721 — that's the reference we need to reach. If that path runs a
DIFFERENT algorithm than sequential branch-by-branch (not just a parallel
implementation), we need to port THAT algorithm, not the sequential loop
we've been porting. Specifically, the BAATARI2 protocol may:
- Collect gains from multiple branches before committing.
- Use a "best-first" acceptance order within an iteration.
- Track `generationofinput[branch]` to avoid re-processing stale branches.
- Do something with `tscorehistory_detail` that our code doesn't.

Effort: 4–8 hours to read the multi-thread path carefully and port its
single-worker behavior to Rust, plus FFI cross-validation to verify we
match at each step.

Alternative (faster, less elegant): instrument `tditeration.c` with
`fprintf(stderr,...)` after each accept to dump `bseq[0]` and its width,
rebuild dvtditr, run one iteration, capture the trace. Compare Rust's
per-branch state against that trace to find the first divergent branch.
Effort: 1–2 hours.

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
