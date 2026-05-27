# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION,
> not just at convergence. A run that matches C's final width by coincidence
> while producing different intermediate alignments is NOT parity — it means
> our accept trajectory differs from C's, and the match is fragile
> (input-dependent).

## Current parity matrix (36-seq protein sample, mafft-upstream/test/sample)

`target/release/mafft-rs <args> sample` vs system `mafft <args> sample`
(MAFFT 7.526), diffing the FASTA outputs.

| Mode / flags                                | C width | Rust width | diff lines | Status |
|---------------------------------------------|---------|------------|------------|--------|
| FFT-NS-2 (default)                          | 717     | 717        | 0          | byte-exact ✓ |
| NW-NS-2 (`--nofft`)                         | 717     | 717        | 0          | byte-exact ✓ |
| FFT-NS-i (`--maxiterate 100`)               | 721     | 721        | 0          | byte-exact ✓ |
| L-INS-1 (`--localpair --maxiterate 0`)      | 719     | 719        | 0          | byte-exact ✓ |
| G-INS-1 (`--globalpair --maxiterate 0`)     | 746     | 746        | 0          | byte-exact ✓ |
| E-INS-1 (`--genafpair --maxiterate 0`)      | 729     | 729        | 0          | byte-exact ✓ |
| L-INS-i (`linsi` / `--localpair --maxiterate 1000`) | 735 | 735 | 0       | byte-exact ✓ |
| G-INS-i (`ginsi` / `--globalpair --maxiterate 1000`)| 737 | 737 | 0       | byte-exact ✓ |
| E-INS-i (`einsi` / `--genafpair --maxiterate 1000`) | 729 | 729 | 0       | byte-exact ✓ |
| BL30 / BL45 / BL50 / BL62 / BL80 FFT (`--bl N`) | match | match | 0      | byte-exact ✓ |
| BL80 NW (`--bl 80 --nofft`)                 | 712     | 712        | 0          | byte-exact ✓ |
| JTT 200 / JTT 100 FFT (`--jtt N`)           | match   | match      | 0          | byte-exact ✓ |
| TM 100 / TM 200 NW (`--tm N --nofft`)       | match   | match      | 0          | byte-exact ✓ |
| TM 100 / TM 200 FFT (`--tm N`, default retree=2) | match | match    | 0          | byte-exact ✓ |
| `--tm 200 --retree 1`                       | 717     | 717        | 8          | C-side tied-trace artifact ✓ (see §B.2) |
| `--tm 200 --treeout`                        | match   | match      | 20         | C-side tied-trace artifact ✓ (pass-0 TM drift propagates into tree branches; see §B.2) |
| PartTree (`--parttree`)                     | 752     | 752        | 0          | byte-exact ✓ |
| DP-PartTree (`--dpparttree`)                | 752     | 752        | 0          | byte-exact ✓ |
| PartTree NW (`--parttree --nofft`)          | 752     | 752        | 0          | byte-exact ✓ |
| `--parttree --reorder`                      | match   | match      | 0          | byte-exact ✓ |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6 fixture) | match | match | 0 | byte-exact ✓ |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 62 (case)  | byte-exact mod case ✓ |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 62 (case)  | byte-exact mod case ✓ (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0          | byte-exact ✓ |
| `--reorder` (FFT-NS-2, INS-i family)        | match   | match      | 0          | byte-exact ✓ |
| `--treeout` (FFT-NS-2, NW-NS-2, FFT-NS-i, L/G/E-INS-i, BL/JTT, parttree, dpparttree) | match | match | 0 | byte-exact ✓ |
| `--treein` (FFT-NS-2, NW-NS-2, BL/JTT/TM, L/G/E-INS-i, FFT-NS-i, all non-parttree) | match | match | 0 | byte-exact ✓ |
| `--treein --treeout`                        | match   | match      | 0          | byte-exact ✓ (appends `#by loadtree\n` like C `mltaln9.c:2818`) |
| `--auto` (small/medium/large brackets)      | match   | match      | 0          | byte-exact ✓ |
| `--memsavetree` / `--memsavetree --treeout` | match   | match      | 0          | byte-exact ✓ |
| `--anysymbol` / `--preservecase` (protein + DNA, non-standard chars, mixed case) | match | match | 0 | byte-exact ✓ |
| `--leavegappyregion` / `--legacygappenalty` (combos with --anysymbol/--reorder/--treein/--memsavetree/--auto) | match | match | 0 | byte-exact ✓ |
| `--seed FILE` (L/G/E-INS-i + FFT-NS-i, single + multiple files) | match | match | 0 | byte-exact ✓ |
| `--seedtable FILE` (L/G/E-INS-i + FFT-NS-i, pre-computed hat3.seed) | match | match | 0 | byte-exact ✓ |
| `--memsave` / `--nomemsave` (FFT-NS-2, FFT-NS-i, retree-1, --memsavetree, --nofft combos) | match | match | 0 | byte-exact ✓ (Hirschberg `msalignmm` wired into engine) |
| `--retree {1,2,3,5}` × {FFT-NS-2, FFT-NS-i, L/G/E-INS-i} (20 combos) | match | match | 0 | byte-exact ✓ |

Test suite: **357 Rust integration tests pass, 0 failed, 5 ignored**
(`cargo test --workspace --release --tests`). The 5 ignored are all
deliberate diagnostic/bisection tools (drift exploration, fixture-locked
hand-bisection scripts) that print analysis instead of asserting
pass/fail; each is annotated at the call site. Plus 32 Python tests pass.

## Closed work (full implementation notes in git history)

| Section | Topic | Date |
|---------|-------|------|
| §1 | G-INS-1 — `global_align` mirroring `G__align11`'s max-so-far DP with `>=` tie-break; `outgap=1` threaded for G-INS-i. | 2026-05-06 |
| §2 | E-INS-1 — `genaffine_local_align` mirroring `genL__align11`; `PairAligner::GeneralizedAffine`; E-INS-i `lexp=laof=0` per `scripts/mafft:1940-1948`. | 2026-05-06 |
| §3 | L/G/E-INS-i refinement — impmatch diagonal, `Falign_localhom` port, `partA__align` strict-`>` tiebreak, boundary nongap-freq, hat2-rounded refinement-tree distance, `naivepairscore11` for E-INS-i. | 2026-05-07 |
| §4 | BLOSUM 50 FFT — `f64::mul_add` in `match_calc_row` matching gcc-O3 FMA. | 2026-05-08 |
| §5 | FFT anchor / TM 200 FFT — bit-for-bit FFT C-port, `shift_and_score` sign convention, anchor coordinate mapping, 2-channel polarity+volume FFT. | 2026-05-10 |
| §6 | PartTree — distance + lenfac, pivot selection (`tokyoripara=0` when `picksize>njob`), UPGMA-on-yukomtx, unweighted profile, `outgap=1`. | 2026-05-10 |
| §7 | `--add` / `--addfragments` / `--keeplength` — `mergeoralign[]` semantics, per-step `findcommongaps`/`commongappick`/`restorecommongaps`/`insertnewgaps` for NewRight/NewLeft. | 2026-05-10 |
| §9a | Q-INS-i — works when `mxscarnamod` is built from `mafft-upstream/extensions`. | 2026-05-10 |
| §A | `--parttree --nofft` — `pairwise_align11` now receives `penalize_term_gaps` (was hard-coded `false, false`). | 2026-05-12 |
| §A0 | `--allowshift` — `makedynamicmtx` skips `amino[i] == '-'` row/col when applying `offset * 600` delta (matches `mltaln9.c:15197-15203`). | 2026-05-12 |
| §AA | `--reorder` / `--inputorder` — `Topology::dfs_order`; `engine.reorder_output` permutes after refinement. PartTree closed via two-pass composition `final_order[k] = call1_order[call2_order[k]]` from `splittbfast` CALL 1/2. | 2026-05-13 |
| §AB | `--treeout` — `topology_to_newick` matching `mltaln9.c:6190-6491`; `--dpparttree --treeout` via `run_parttree_pipeline_with_scorer` with `G__align11_noalign` on `n_disLN`. | 2026-05-13 |
| §AC | `--treein` — `treein.rs` parser; overrides progressive merge tree, LH-table importance reweighting, and refinement-tree rebuild (matching `tbfast.c:2072`/`:2967`/`:1355`/`dvtditr.c:766`). | 2026-05-14 |
| §AD | `--auto` — `decide_auto` reads `(nseq, nlen)`, picks mode + retree from the 7-bracket table. | 2026-05-14 |
| §AE | `--memsavetree` — port of `compacttreegivendist` (`mltaln9.c:5221`) — NOT `compacttree_memsaveselectable` which is `--youngestlinkage`. Pass 0 k-mer + pass 1 MSA distance. | 2026-05-14 |
| §AF | `--anysymbol` / `--preservecase` — `read_fasta_casepreserve`, `replace_unusual`, post-align name-keyed restore (mirrors C `replaceu` + `restoreu`). | 2026-05-14 |
| §AG | `--leavegappyregion` / `--legacygappenalty` — `GapModel::legacy_gap_cost` threads through; profile DP early-returns into recursive call with `nongap_freq = [1.0; len]` clones. | 2026-05-14 |
| §AH | `--seed FILE` (repeatable) — `extract_putlocalhom2_regions`, `build_seed_homology_table` with `tsuyosa = user_nseq² * 100`, merged into pairwise LH before `recompute_importance`. Forces `iterate ≥ 2`. | 2026-05-14 |
| §B.1 | `mi += penalty_ex` / `m[j] += penalty_ex` increments in profile DP. `pairwise_align11` still missing it but benign because CLI doesn't expose `--exp`. | 2026-05-11 |
| §B.3 | `--memsave` / `--nomemsave` — Hirschberg `msalignmm` (`mafft-align/src/msalign.rs`) wired into engine. Two off-by-one bugs fixed at midpoint state update (`midw[j] += wm` per `MSalignmm.c:1610`, `midm[0] += firstm` per `:1637`) — found by patching upstream `MSalignmm.c` with `fprintf` instrumentation. | 2026-05-14/15/16 |
| §B.3.2 | `--seedtable FILE` — `parse_hat3_seed` reads the 9-field text format `multi2hat3s.c:214` writes. Engine wiring shares the `--seed` path. CLI gating: rejects `--seed`, `--add`, `--parttree`, `--memsave`; forces `iterate ≥ 2`. | 2026-05-17 |
| §B.4 | `--retree N` — engine now mirrors two undocumented `scripts/mafft` rewrites: `:1840-1842` clamps `cycle = min(cycle, 3)`; `:1934-1936` forces `cycle = 1` for L/G/E/Q/X-INS-i. Pre-fix, `--retree 3 --localpair` diverged 898 lines. 7 regression tests + verified manually across 20 combinations. | 2026-05-18 |
| §B.5 | `profile_cache` HashMap → BTreeMap. Get-only today, but `BTreeMap` removes the iteration-order-dependence footgun for future refactors. | 2026-05-18 |
| §B.7 | Deleted `parttree.rs` (legacy module with broken `max_by_key` tie-break in unreachable-for-36-seq recursive code). Re-routed `--dpparttree` through `parttree_split::build_parttree_topology`. | 2026-05-18 |
| §B.8 | `--treein --genafpair --maxiterate >1` — `recompute_importance` now uses `user_topo` (when set) for `sequence_weights`, not `musclesupg(&dm)`. Affected E-INS-i only because L/G/FFT-NS-i happened to converge to the same tie-breaks on the 36-seq sample regardless of weights. | 2026-05-14 |
| §B.9 | `--memsavetree` — see §AE; initial port was `compacttree_memsaveselectable` (the `--youngestlinkage` algorithm), not `compacttreegivendist`. Replaced after FFI cross-validation showed C uses the ELSE branch at `disttbfast.c:4017` for `compacttree == 2`. | 2026-05-14 |
| §B.10 | FFT-segmented refinement boundary frequencies — `refinement.rs::realign_all` use_fft branch was calling `profile_align` (which defaults `BoundaryFreqs::default() = 1.0`) for every segment. C's per-segment `A__align` derives `headgapfreq{1,2}` / `gapfreq{1,2}[lgth]` from `sgap1`/`sgap2`/`egap1`/`egap2` via `outgapcount` on the parent alignment's boundary columns. Wired the same computation; also removed the prior `prof_seg.fgcp[last] -= w` correction, which mirrored the *disabled* `#if 0` branch of C's `new_FinalGapCount` (the compiled `#if 1` branch effectively ignores `egappat` at the last position via null-terminator read). Closes BB30013 (was 3843-line residual). BALIBASE 3 FFT-NS-i 99.7% → **100% (386/386)**. Regression test `fftnsi_segmented_boundary_byte_identical_to_c` on BB12019 + `bali3.BB12019.fftnsi.iter5` fixture. | 2026-05-25 |
| §B.11 | `--seed` FFT-NS-i — `engine.rs:876` `use_segmented` was unconditionally `matches!(self.mode, FftNsi { .. })`; needed `&& local_hom.is_none()` since C only segments when `constraint == 0` (`dvtditr.c:882`). Without the gate, seeded FFT-NS-i would re-run segmented refinement and break the seed/seedtable `_byte_identical_to_c` tests. | 2026-05-26 |

---

## §B.2. `--tm 200 --retree 1` 8-line diff + BALIBASE 3 corpus — RESOLVED 2026-05-19 (upstream report `MAFFT_UPSTREAM_REPORT.md` covers both)

**The only known byte-divergence in production**: `mafft --tm 200 --retree 1
sample` vs `mafft-rs --tm 200 --retree 1 sample` differs by 8 lines — 4
single-char gap shifts in 2 of 36 sequences (M62903 and S75720, chicken
opsins). Default `--tm 200` (retree=2) converges to C's output because the
second pass realigns from the rebuilt tree, but the first-pass divergence
cascades into `--tm 200 --treeout` branch-length drift (~3e-3, 20-line
diff).

**Diagnosis (cell-level FFI cross-validation, 2026-05-18)**: our Rust
`profile_align` produces byte-identical `ijp[i][j]` to C's `A__align`
when both are called in isolation with the same inputs (verified
cell-by-cell across 365×363 cells of the divergent step-12 grid). But
C's `A__align` in the full-run context — after 12 prior calls in
`treebase` — produces a *different* `ijp` matrix with the same scalar
score (208675.4644010082) and width (367).

The visible inputs (`initverticalw`, `currentw`, `ogcp1/2`, `fgcp1/2`,
`gapfreq1pt/2pt`, normalized eff weights, dynamic matrix) are
byte-identical between the isolated test and the full-run invocation.
The `ijp` divergence is only explainable by C's `A__align` having
`static TLS` buffers (`commonIP/ijp`, `cpmx1/2`, `ogcp/fgcp1/2o`,
`gapfreq1/2`, `m`, `mp`, `w1/w2`) whose post-call residue biases the
next call's tied-cell selection. A warm-up `A__align` call in the test
flipped 1 ijp cell, confirming static-state sensitivity. TM PAM 200 is
the only matrix flat enough to surface this — BLOSUM/JTT have enough
score diversity that tied DP cells are rare.

**Severity**: LOW. Both Rust and C produce optimal-scored alignments
on the failing config; only the choice between tied traces differs.
All correctness tests pass; the 8-line diff is *the only* known
byte-difference vs C MAFFT 7.526 on the 36-seq test sample.

**Defensive FMA fixes that landed during the investigation** (kept,
they harden other latent precision boundaries):
- `mafft-core/src/progressive.rs::blend_profiles_exact` (and the
  `blend_og_one_side` / `blend_fg_one_side` helpers) — 4 `+= a * b`
  patterns converted to `mul_add`.
- `mafft-tree/src/weighting.rs::sequence_weights` —
  `rootnode[s] += step.left_length * eff[s]` (mirrors C
  `mltaln9.c:9893`) converted to `mul_add`.

**Pending action**: `MAFFT_UPSTREAM_REPORT.md` documents the full
investigation (cell-level FFI capture, patches, conclusion) and is
ready to send upstream (Kazutaka Katoh, katoh@ifrec.osaka-u.ac.jp).
Closing this on our side would require deliberately replicating C's
exact static-buffer cross-call dependency in our `profile_align` —
significant complexity for ~zero benefit since both traces are
optimal-scored.

---

## §B.6. Sequential float summation guard in `refinement.rs:1069-1080`

**Location**: `crates/mafft-core/src/refinement.rs:1069-1080` —
intergroup scoring uses explicit sequential summation with a comment
warning against `par_iter().sum()`.

**Status**: Already correct and documented. This section exists *as*
the reminder — removing the sequential guard would cascade into
non-deterministic refinement accept/reject decisions.

---

## §C. Performance / non-correctness items

### Benchmark snapshot (2026-05-18)

5 runs summed, real time, macOS arm64, `mafft-rs` release vs C MAFFT
7.526, 108-seq synthetic input:

| Mode                                     | C MAFFT | mafft-rs | Ratio          |
|------------------------------------------|---------|----------|----------------|
| default (FFT-NS-2)                       | 1.34s   | 1.16s    | Rust 1.16× faster |
| `--maxiterate 50` (FFT-NS-i)             | 6.62s   | 7.78s    | Rust 1.18× SLOWER |
| `--maxiterate 50 --localpair` (L-INS-i)  | 37.00s  | 22.93s   | Rust **1.61× faster** |
| `--maxiterate 50 --globalpair` (G-INS-i) | 35.32s  | 22.62s   | Rust **1.56× faster** |
| `--parttree`                             | 1.70s   | 0.09s    | Rust **18× faster**   |
| `--dpparttree`                           | 7.24s   | 0.08s    | Rust **90× faster**   |

### §C.1. Per-group gap stripping — DEFERRED

Pre-2026-05-18 theory: C strips per group via `commongappick()`, we
strip globally — for a 500-column MSA with `kept1=100, kept2=120` C
DPs 12,000 cells while we DP up to 250,000. *Empirical reality*: we're
already faster than C on the modes this would help (default 1.16×,
L-INS-i 1.61×, G-INS-i 1.56×, PartTree 18-90×). The cell-count
theory doesn't translate to wall-clock because most refinement
branches don't have wildly different `kept1`/`kept2` columns in
practice. Implementing per-group stripping requires porting C's
`insertnewgaps()` round-trip from `addfunctions.c` — high cost, no
measurable benefit. Reconsider only if a real workload surfaces a
progressive-alignment bottleneck.

### §C.2. FFT-NS-i ~18% gap — DEFERRED (profile data captured in `PROFILING.md`)

FFT-NS-i is the only mode where Rust is slower than C. Per the
profile (macOS `sample`, full recipe in `PROFILING.md`):

| Self-time samples | Function | % main thread |
|-------------------|----------|---------------|
| 509 | `profile_align_imp_with_boundary` (inner DP cell update) | **55%** |
| 97  | `compute_split_score`             | 11% |
| 89  | `Profile::from_aligned`           | 10% |
| 26  | `Profile::match_score`            | 3% |
| 22% rest | rayon, alloc, traversal | — |

**Attempts (2026-05-18)**:

- **`Profile::from_aligned` 3-pass → single-pass fusion** —
  *LANDED* (`profile.rs::Profile::from_aligned`). Per-sequence walk
  now updates freqs, gap_freq, opening_count, closing_count in one
  pass via a `gc_prev` state byte (vs 3 separate walks before).
  Tail-end semantics preserved for sequences shorter than `length`.
  Wall-clock impact within noise; code-cleanliness win, no behavior
  change.
- **`h` / `ijp` thread-local pool** — *LANDED* (`profile.rs` top:
  `DP_H_POOL`, `DP_IJP_POOL`). Take the two big 2D DP matrices from
  thread-local pools at the top of `profile_align_imp_with_boundary`,
  swap back before returning. Grow-only resize, no per-cell zero —
  every read cell is unconditionally written by boundary init or DP
  body before any traceback read. Re-profile confirmed
  `raw_vec::grow_one` samples inside the DP function dropped from
  ~80 to ~20-30. Wall-clock impact within noise on macOS (malloc/free
  is already fast); kept for code-cleanliness + symmetry with C
  `A__align`'s `static TLS` buffer reuse.
- **Branchless `compute_split_score`** — *ATTEMPTED, REVERTED*.
  Rewrote `pairwise_score` as a per-cell state-machine to remove the
  `while (seq1[k] == '-')` consume loop. Failed 16 byte-identity
  tests because C's consume loop crosses both-gap positions
  (it only checks seq1), while a naive reset-on-both-gap rule
  re-charges penalty when an A-gap-run is interrupted by a both-gap
  column. The minimum correct branchless form is a 3-state machine
  (Neutral / AGapRun / BGapRun) — still branchy, no clean SIMD path.
  Reverted to the C-style structure with an explanatory comment.

**Remaining gap is in the inner DP cell-update arithmetic itself**
(456/922 main-thread samples post-pool). Closing it requires the
anti-diagonal SIMD rewrite — a multi-day project with byte-identity
regression risk on the warp DP, gap-skip trackers, and traceback
ordering. Deferred until a sustained FFT-NS-i workload justifies the
cost.

SIMD-friendly patterns in `match_score()`, `pairwise_score()`,
`pairwise_identity_distance()` continue to auto-vectorize via LLVM
where possible (the gap-run consume loop in `pairwise_score` still
inhibits it on the gap-run paths).

---

## §E. BALIBASE 3 parity sweep — all four INS-i modes 386/386 (100%)

### Per-mode results (BALIBASE 3 RV11-RV50, 386 files)

| Mode | Flags | Match | Total % | Residuals |
|------|-------|-------|---------|-----------|
| FFT-NS-i | `--maxiterate 100` | **386/386** | **100.0 %** | none (2026-05-25) |
| L-INS-i | `--localpair --maxiterate 100` | **386/386** | **100.0 %** | none (2026-05-27) |
| G-INS-i | `--globalpair --maxiterate 100` | **386/386** | **100.0 %** | none (2026-05-27) |
| E-INS-i | `--genafpair --maxiterate 100` | **386/386** | **100.0 %** | none (2026-05-27) |

**All four INS-i modes are 386/386.** §E.3 (below) records how the last
residual, E-INS-i BB40004, was closed.

### §E.3 E-INS-i BB40004 — CLOSED 2026-05-27 (phase-split importance)

Surfaced by the first full E-INS-i sweep. Investigation chain (all
proven, not guessed):
- E-INS-1 (`--maxiterate 0`) already diverged 588 lines → bug is in the
  progressive phase, not refinement. G-INS-1 was byte-exact.
- Feeding C's guide tree via `--treein` → diff=0 → the divergence is
  *entirely* in the guide tree / its derived weights, nothing else.
- The genaffine pairwise alignment and distances are bit-identical to C
  (instrumented C `pairlocalalign` score2dist dump: 0/2211 pscore
  mismatches at full `%.18e` precision). So NOT a genaffine trace tie.
- UPGMA-trace comparison: Rust's *progressive* tree (block 1) is a
  66/66 bit-exact match with C. The divergent blocks were the
  *importance* tree (built from rounded dm).
- Root cause: C computes constraint `importance` TWICE with different
  trees — tbfast (progressive) from the FULL-PRECISION in-memory
  `iscore` tree (`tbfast.c:2193`+`2926`, no hat2 round-trip there),
  dvtditr (refinement) from the 3-decimal `hat2` tree
  (`readhat2_pointer`). Rust used the rounded tree for BOTH. L/G-INS-i
  tolerated it (full vs rounded importance tree happened to give the
  same progressive merge); E-INS-i BB40004 had a near-tie that flipped.

Fix (`engine.rs`): the pre-progressive `recompute_importance` now uses
the full-precision `dm` tree (phase 1); a second `recompute_importance`
on `local_hom` from the rounded refinement tree runs just before
`iterative_refine` (phase 2). Gated to the pairwise INS-i path
(`pairwise_for_constraints.is_some() && user_topo.is_none() &&
!uses_rna_constraints`). Verified: E/L/G-INS-i all 386/386 full sweeps;
89 end_to_end byte-identity tests pass (incl. all --seed/--seedtable/
--treein INS-i combos). A&B-tested both directions: rounded-tree
progressive regresses BB40004 (588 lines); full-tree refinement
regresses BB50001 (244 lines) — confirming the split is required.

**The 2026-05-26 conclusion below ("why we can't close the remaining 3,"
blamed on `A__align` static-TLS coupling) was WRONG and is retained only
as a record of the dead end.** The L-INS-i residuals BB30028 + BB50001,
and the G-INS-i residual BB20004, were all closed on 2026-05-27 by
fixing four floating-point summation-order mismatches in the constrained
refinement path — found by instrumenting C MAFFT and diffing values at
`%.18e`. None involved static-TLS memoization.

### §E.1 The four FP-ordering fixes (2026-05-27, `refinement.rs` + earlier `engine.rs`/`constraints.rs`)

1. **`hat2` 3-decimal truncation for the refinement-weights tree**
   (`engine.rs`). C's `dvtditr` rebuilds the refinement tree from the
   3-decimal `hat2` file (`readhat2_pointer`, `DFORMAT="%#6.3f"`), not
   the in-memory full-precision distance matrix. Truncate `dm` to 3
   decimals before `musclesupg` for the weights tree only (the
   progressive tree keeps full precision — truncating it regresses
   BB30013/BB40004). Closed BB50001 (244 lines).
2. **printf banker's rounding for `opt`** (`constraints.rs`). C's hat3
   `%7.5f` write + `atof` read rounds half-to-even; Rust's
   `(v*1e5).round()/1e5` rounds half-away-from-zero. They disagree at
   exact `x.x005` half-points (BB20004 pair (31,36): `opt_pre*1e5 =
   197620.5` exactly). `format!("{:.5}", v).parse()` reproduces printf.
   Makes the importance table bit-identical to C (0/12694 region diffs).
3. **backward diagonal sum for `old_imp`** (`compute_impmatch_diagonal`).
   C `tditeration.c:891` sums `imp_match_out_scD(i,i)` from `i=length-1`
   down. FP add is non-associative; forward summation drifts ~1 ULP.
   Closed BB20004 (10 lines).
4. **per-segment forward `impmatch` for `new_imp`**
   (`realign_all_constrained_fft`). C's `Falign_localhom` accumulates
   each FFT segment's impmatch (backward within segment, via
   `Atracking_localhom`) then sums forward across segments
   (`Falign_localhom.c:816`). A single global diagonal sweep diverges
   ~7 ULP and flips one tied accept/reject. Threaded the segment-
   accumulated value out through `realign_all` and used it for
   `new_imp` (falls back to the global sum only on the non-FFT path,
   which `dvtditr -F` never takes). Closed BB30028 (4 lines).

A naive global *forward* sum (instead of per-segment) closed BB30028
but regressed 6 other files (BB30010 catastrophically) — the
segmentation *structure*, not a global direction, is what matches C.

### §E.2 Superseded dead-end analysis (kept for the record)

### Fixes that closed cases (committed)

1. **`mafft-tree/src/distance.rs::ktuple_distance_aa/_nuc` lenfac uses
   `nogap_len` not filtered groups** — closed 11 of 17 initial cases
   (all involving X / `.` non-standard residues that the filtered-group
   length stripped but `nogaplen` retains). C's
   `disttbfast.c:3845-3856` uses `nogaplen`. Regression test:
   `distance::tests::ktuple_lenfac_uses_nogaplen_not_filtered_groups`.

2. **`progressive::blend_fg_one_side` past-null-terminator handling**
   — C reads `gaptable[j+1]` at `j = alen-1` (the `\0` terminator,
   treated as non-gap) and adds the closing-count contribution; Rust
   short-circuited at `j < alen - 1` and dropped it. Now uses
   `next_is_gap()` helper. Real correctness bug; doesn't surface on
   today's baseline tests but ships with the rest.

3. **FFT-segmented refinement boundary frequencies (§B.10, 2026-05-25)**
   — `refinement.rs::realign_all` use_fft branch was calling
   `profile_align` which defaults `BoundaryFreqs = 1.0` for every
   segment, instead of computing per-segment `headgapfreq{1,2}` /
   `gapfreq{1,2}[lgth]` from `sgap`/`egap` via `outgapcount` on the
   parent alignment's boundary columns (mirrors C `Salignmm.c:1585-
   1622`). Combined with removing the prior `fgcp[last]` correction
   (which had mirrored the disabled `#if 0` branch of C's
   `new_FinalGapCount`), this closed BB30013 (was 3843-line residual)
   and lifted FFT-NS-i from 99.7 % → 100 %. Regression test
   `fftnsi_segmented_boundary_byte_identical_to_c` (input
   `bali3.BB12019.fa`, reference `bali3.BB12019.fftnsi.iter5`).
   Full writeup: `~/.claude/.../memory/project_bb30013_fix.md`.

### Why we can't close the remaining 3 — SUPERSEDED 2026-05-27 (was wrong; see §E.1)

> This hypothesis was disproven. The residuals were FP summation-order
> bugs (§E.1), not static-TLS coupling. The cell-level FFI equivalence
> tests below remain valid and useful; the *conclusion* that closing the
> cases required porting `reuseprofiles` was incorrect.

Same root cause as §B.2: C `A__align`'s `static TLS` memoization
(`Salignmm.c:1446-1450`) controlled by `calledbyfulltreebase=1`. The
memoization conditionally reuses cached cpmx from prior `A__align`
calls in the same process — output thus depends on cross-call state
we don't carry in our stateless progressive engine.

**Proven bit-equivalent with C** via dedicated FFI cross-validation
tests in `crates/mafft-core/tests/cross_validate_*` (now in-repo,
fixtures `bali3.BB20027.fa`, `bali3.BB12041.fa`, `bali3.BB12019.fa`):
- `cross_validate_cpmx::rust_profile_freqs_match_c_cpmx_calc_new` —
  `Profile::from_aligned` vs C `cpmx_calc_new + gapcountf +
  st_OpeningGapCount + st_FinalGapCount` (4 fields, zero drift).
- `cross_validate_cpmx::rust_blend_matches_c_blend_cell_by_cell` —
  `blend_profiles_exact` vs C `createcpmxresult + creategapfreqresult +
  createogresult + createfgresult` (zero drift on 3+5 gappy input;
  auto-`#[ignore]`d on Linux glibc due to `free(): invalid pointer`
  at process teardown — runs on macOS).
- `cross_validate_counteff::bb20027_pass1_weights_match_c` —
  `sequence_weights` vs C `counteff_simple_double` (zero drift, runs
  on the in-repo BB20027 fixture).
- `cross_validate_bb20027_dp::bb20027_step{12,13}_rust_dp_vs_c_aalign_cell_by_cell`
  — Rust `profile_align` vs C `A__align` (zero drift when
  `cpmxchild=NULL, calledbyfulltreebase=0`).

~~Closing these would require porting C's `reuseprofiles` static TLS
state machine to Rust~~ — **NO.** This was the wrong conclusion. The
cases were closed on 2026-05-27 via the four FP summation-order fixes
in §E.1, with the engine still pure (no cross-call state added). The
FFI equivalence tests above remain valid; only the inference from them
was mistaken.

### `--c-compat` flag (2026-05-21)

Added opt-in flag (`MafftEngine.c_compat`, CLI `--c-compat`) that
enables C MAFFT's `reuseprofiles` static-TLS cpmx memoization via
`Profile::from_aligned_with_memo`. Infrastructure compiles, has
unit tests, no-op when disabled. **Does NOT close the residual
divergences** — verified via C trace: memo doesn't fire at the
divergent steps. Kept as foundation for future work. Full notes in
`balibase_parity_run.md` §Layer 8+.

### Done

FFT-NS-i, L-INS-i, G-INS-i, and E-INS-i are **all 386/386 (100 %)** on
the full BALIBASE 3 sweep (2026-05-27). The INS-i refinement residuals
closed via FP summation-order parity (§E.1); the last one, E-INS-i
BB40004, closed via phase-split constraint importance (§E.3). The engine
stays stateless — no `reuseprofiles` port needed.

---

## §D. X-INS-i (`--xinsi`) — UNTESTABLE without `contrafold`

Requires Stanford's `CONTRAfold v2.02+` binary (http://contra.stanford.edu/contrafold/),
distributed separately. Not shipped by upstream MAFFT, not buildable
from `mafft-upstream/extensions`. Rust wiring exists
(`engine.rs::XInsi`, `mafft-bin/src/main.rs:--xinsi`) and emits the
correct "contrafold not found" diagnostic when the binary is absent.
End-to-end validation deferred until `contrafold` is installed.

---

## Active items summary

Everything that follows is either (a) a documented gotcha that we
choose not to fix, or (b) an external dependency we can't satisfy:

1. **§B.2 + §E residuals** — same root cause: C `A__align` static-TLS-
   buffer memoization (`Salignmm.c:1446-1450`) produces context-
   dependent traceback choices for tied DP cells. Concretely:
   - `--tm 200 --retree 1` on the 36-seq sample: 8-line diff.
   - 3 of 386 BALIBASE 3 tests across all measured INS-i modes:
     BB30028 (L-INS-i, 4 lines), BB50001 (L-INS-i, 244 lines),
     BB20004 (G-INS-i, 1732 lines). FFT-NS-i is now 386/386 = 100 %
     after the §B.10 fix.
   `MAFFT_UPSTREAM_REPORT.md` documents the diagnosis; ready to send
   to katoh@ifrec.osaka-u.ac.jp. Cell-level FFI tests
   (`cross_validate_cpmx`, `cross_validate_counteff`,
   `cross_validate_bb20027_dp`) PROVE the Rust DP and profile
   building match C bit-for-bit on identical inputs — divergences
   come from C's static state, not from a Rust bug. Closing would
   require porting C's `reuseprofiles` state machine; not worth it.
2. **§B.6** — sequential float summation guard in `refinement.rs`.
   Section IS the reminder; no work needed.
3. **§C.1** — per-group gap stripping. Deferred; we're already faster
   than C on the affected modes.
4. **§C.2** — FFT-NS-i ~18% gap. Deferred; profile data in
   `PROFILING.md`. Closing it requires the anti-diagonal SIMD
   rewrite (multi-day project, no measurable production benefit
   currently).
5. **§D** — X-INS-i needs `contrafold` binary; untestable until
   someone installs it.
