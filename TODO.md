# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION,
> not just at convergence. A run that matches C's final width by coincidence
> while producing different intermediate alignments is NOT parity — it means
> our accept trajectory differs from C's, and the match is fragile
> (input-dependent).

## Current parity matrix (36-seq protein sample, mafft-upstream/test/sample)

Verified 2026-05-12 by running `target/release/mafft-rs <args> sample` against
the system `mafft <args> sample` (MAFFT 7.526) and diffing the FASTA outputs.
Widths are taken as `length(seq[1])` (the canonical MSA width); the system
mafft and our binary agree on every byte for every ✓ row.

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
| BL30 FFT (`--bl 30`)                        | 773     | 773        | 0          | byte-exact ✓ |
| BL45 FFT (`--bl 45`)                        | 729     | 729        | 0          | byte-exact ✓ |
| BL50 FFT (`--bl 50`)                        | 712     | 712        | 0          | byte-exact ✓ |
| BL62 FFT (`--bl 62`)                        | 717     | 717        | 0          | byte-exact ✓ |
| BL80 FFT (`--bl 80`)                        | 712     | 712        | 0          | byte-exact ✓ |
| BL80 NW (`--bl 80 --nofft`)                 | 712     | 712        | 0          | byte-exact ✓ |
| JTT 200 FFT (`--jtt 200`)                   | 729     | 729        | 0          | byte-exact ✓ |
| JTT 100 FFT (`--jtt 100`)                   | 732     | 732        | 0          | byte-exact ✓ |
| TM 100 NW (`--tm 100 --nofft`)              | 767     | 767        | 0          | byte-exact ✓ |
| TM 200 NW (`--tm 200 --nofft`)              | 765     | 765        | 0          | byte-exact ✓ |
| TM 100 FFT (`--tm 100`)                     | 767     | 767        | 0          | byte-exact ✓ |
| TM 200 FFT (`--tm 200`)                     | 767     | 767        | 0          | byte-exact ✓ |
| PartTree (`--parttree`)                     | 752     | 752        | 0          | byte-exact ✓ |
| DP-PartTree (`--dpparttree`)                | 752     | 752        | 0          | byte-exact ✓ |
| PartTree NW (`--parttree --nofft`)          | 752     | 752        | 0          | byte-exact ✓ (closed 2026-05-12) |
| `--add` (30+6 fixture)                      | 741     | 741        | 0          | byte-exact ✓ |
| `--add --nofft` (30+6 fixture)              | 741     | 741        | 0          | byte-exact ✓ |
| `--add --keeplength` (30+6 fixture)         | 595     | 595        | 0          | byte-exact ✓ |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 62 (case)  | byte-exact mod case ✓ |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 62 (case)  | byte-exact mod case ✓ (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0          | byte-exact ✓ (closed 2026-05-12) |
| `--reorder` (FFT-NS-2, INS-i family)        | match   | match      | 0          | byte-exact ✓ (closed 2026-05-13) |

Test suite as of 2026-05-13: **271 Rust tests pass, 0 failed, 0 ignored**
(`cargo test --workspace --exclude pymafft --release`). Plus 32 Python tests
pass. Every mainstream mode in the matrix above is byte-identical to C
MAFFT 7.526.

Resolved sections (full implementation notes in git history):
- §1 G-INS-1 (2026-05-06) — `global_align` ported to mirror `G__align11`'s
  max-so-far DP with `>=` tie-break; `outgap = 1` threaded for G-INS-i.
- §2 E-INS-1 (2026-05-06) — `genaffine_local_align` ported to mirror
  `genL__align11`; `PairAligner::GeneralizedAffine` added; E-INS-i
  parameter overrides (lexp=laof=0) mirrored from `scripts/mafft:1940-1948`.
- §3 L/G/E-INS-i refinement (2026-05-07) — six combined fixes: impmatch
  diagonal, Falign_localhom port, partA__align strict-`>` tiebreak,
  boundary nongap-freq, hat2-rounded distance for refinement tree,
  naivepairscore11 for E-INS-i distance.
- §4 BLOSUM 50 FFT (2026-05-08) — `f64::mul_add` for FMA-equivalent
  rounding in `match_calc_row`, matching gcc-O3 FMA output bit-for-bit.
- §5 FFT anchor / TM 200 FFT (2026-05-10) — bit-for-bit FFT C-port,
  `shift_and_score` sign convention, anchor coordinate mapping,
  `align_with_anchors` as segment-boundary DP, 2-channel polarity+volume
  FFT for protein.
- §6 PartTree (2026-05-10) — distance + lenfac, pivot selection
  (tokyoripara=0 when picksize>njob), UPGMA-on-yukomtx, unweighted
  profile mode, `outgap = 1` for parttree.
- §7 `--add` / `--addfragments` / `--keeplength` (2026-05-10) —
  `mergeoralign[]` semantics (iter().all not any), per-step
  findcommongaps/commongappick/restorecommongaps/insertnewgaps for
  NewRight/NewLeft merge types.
- §9a Q-INS-i (2026-05-10) — works when `mxscarnamod` is built from
  `mafft-upstream/extensions`.

---

## §AA. `--reorder` — RESOLVED 2026-05-13 (byte-identical to C, except PartTree)

**Mode**: `mafft --reorder sample` now produces Rust output byte-identical
to C MAFFT 7.526 across all non-PartTree modes (FFT-NS-2, FFT-NS-i,
NW-NS-2, L/G/E-INS-1, L/G/E-INS-i). Verified 2026-05-13 on the 36-seq
sample.

### Implementation

1. **`crates/mafft-tree/src/topology.rs::Topology::dfs_order`** — returns
   leaves in tree-DFS order. Mirrors C's `topolorderz` (`mltaln9.c:1928`):
   each merge step's `left` / `right` vectors already accumulate subtree
   leaves in DFS order, so the root step's `left ++ right` is the full
   ordering.
2. **`crates/mafft-core/src/engine.rs`** — new `reorder_output` field +
   `with_reorder()` setter. Captures the final progressive guide tree
   (`final_progressive_topo`) on each retree pass; after refinement,
   permutes `msa.sequences` and `msa.names` via `topo.dfs_order()`.
3. **`crates/mafft-bin/src/main.rs`** — adds `--reorder` and
   `--inputorder` CLI flags (mutually exclusive). `--inputorder` is the
   default and is a CLI-only flag (no engine effect).

### PartTree caveat (not addressed)

PartTree skipped via `if self.reorder_output && !use_parttree`. Two
reasons:
1. C's `splittbfast` emits the order from its initial recursive
   partition tree (`splittbfast.c:3049`), not from any global UPGMA
   DFS — partition-discovery order, not tree-DFS order.
2. Our `parttree_split::assemble_topology` sorts leaves at each step
   (`parttree_split.rs:161`), so even the topology DFS would not
   reconstruct the discovery order.

Closing this would require either tracking discovery order separately
or reordering inside `assemble_topology`. Listed as residual gap.

---

## §A. `--parttree --nofft` — RESOLVED 2026-05-12 (byte-identical to C)

**Mode**: `mafft --parttree --nofft sample` now produces Rust width
**752** matching C width **752** byte-for-byte (0 diff lines).

### Root cause

`crates/mafft-core/src/progressive.rs::merge_step_cached` called
`pairwise_align11` for the 1-vs-1 no-constraint case with `head_gap` and
`tail_gap` **hardcoded to `false, false`**, ignoring the
`penalize_term_gaps` flag set by `outgap = 1` for `--parttree`
(`scripts/mafft:2655`). C MAFFT terminal-gap-penalizes 1-vs-1 merges
under `--parttree`, so every PartTree NW pair-merge diverged. FFT path
was unaffected because it routes through a different 1-vs-1 helper.

### Fix

`crates/mafft-core/src/progressive.rs::merge_step_cached` — pass
`penalize_term_gaps, penalize_term_gaps` instead of `false, false`:

```rust
pairwise_align11(
    &aligned[group1[0]], &aligned[group2[0]],
    &scoring.consweight_matrix, &scoring.amino_map,
    scoring.gap.open as f64, penalize_term_gaps, penalize_term_gaps,
)
```

All other PartTree variants (`--parttree` FFT, `--dpparttree`) were
already byte-identical; this fix touches only the NW-no-constraint path.

---

## §A0. `--allowshift` — RESOLVED 2026-05-12 (byte-identical to C)

**Mode**: `mafft --allowshift --globalpair --maxiterate 0 sample` now
produces Rust width **1029** matching C width **1029** byte-for-byte
(0 diff lines).

### Root cause (final)

C's `makedynamicmtx` (`mltaln9.c:15197-15203`) applies the `offset * 600`
delta to the substitution matrix BUT SKIPS the row and column where
`amino[i] == '-'` (gap character — index 24 in the protein alphabet).
Rust's `make_dynamic_matrix` (`progressive.rs`) and the per-pair
`dyn_matrix` construction (`constraints.rs`) were applying delta to
EVERY cell, including the '-' row/col.

For raw amino-acid input the DP never directly indexes the '-' row of
the substitution matrix via `score_at(c1, c2)`. But C's static
`amino_dynamicmtx` is char-indexed and the unshifted '-' row's values
feed into the boundary handling of `match_calc_mtx` at the last row
(i = lgth1) when C reads `seq1[0][lgth1] = '\0'` and looks up
`amino_dynamicmtx['\0'][...]`. Without the '-' skip, our shifted matrix
gave a different boundary score from C's, which propagated through the
warp DP and the per-step / per-pair re-align paths, producing the
TRGP-vs-TRG--P alignment shift in pair (U22180, M62903).

### Final fix

1. **`progressive.rs::make_dynamic_matrix`** — added `gap_idx` parameter
   and skip `if i == gap_idx || j == gap_idx { v }` matching C's
   `amino[i] == '-'` check. The caller (`progressive_align_full`)
   computes `gap_idx = scoring.amino_map[b'-' as usize]` once and threads
   it in.
2. **`constraints.rs::build_homology_table_with_unalign`** — per-pair
   dynamic matrix construction in the re-align branch now applies the
   same '-' skip.
3. **Engine-side script-param zeroing** (`engine.rs:304-318`) — restored
   `lexp = laof = 0` when `unalign_active`, mirroring `scripts/mafft:1469-1473`.
   This was previously left off because of a downstream bug; with that
   downstream bug now fixed, this is the correct behavior.

### Trajectory (closed)

Rust progression on n=36 with `--allowshift --globalpair --maxiterate 0`:
- 746 (no effect, only CLI flag wired)
- 809 (per-step matrix scaling, int matrices)
- 957 (+ per-pair re-align, int matrices)
- 993 (after f64 DP migration)
- 987 (after warp DP port)
- 1018 (after `§B.1` `mi += penalty_ex` fix)
- **1029 = C byte-identical** (after `'-'` row/col skip in `make_dynamic_matrix`)

---

## §B. Latent divergence candidates (potential future surfacing)

These are byte-identical on all 36-seq-sample modes but may surface for
other inputs / mode combinations. Investigate when adding a new mode or
when a divergence shows up in the wild.

### §B.1. ~~Missing `mi += penalty_ex` / `mj[j] += penalty_ex` in profile DP~~ — RESOLVED 2026-05-11

**C reference**: `mafft-upstream/core/Salignmm.c:1718, 1727, 1933, 1953`.

C does four `fpenalty_ex` accumulations that Rust was missing:
1. `initverticalw[i] += fpenalty_ex * i` (line 1718, head_gap init).
2. `currentw[j] += fpenalty_ex * j` (line 1727, head_gap init).
3. `mi += fpenalty_ex` (line 1933, in main DP loop after mi update).
4. `if (j < lgth2) m[j] += fpenalty_ex` (line 1953, after mj update).

**Fix landed** in `crates/mafft-align/src/profile.rs::profile_align_imp_with_boundary`:
- Lines 535, 567: init increments.
- Lines 671, 688-690: per-cell increments in main DP loop.

**Guarded by** `profile_align_imp_nonzero_penalty_ex_matches_c_a__align`
in `crates/mafft-core/tests/cross_validate_profile_align.rs`, which runs
profile_align_imp with `penalty_ex = -100` against C's `A__align` with
the same and asserts byte-identical width AND score. Before fix: Rust
score 107421 vs C 104721 (Δ=2700). After fix: byte-identical.

**Why it was latent**: For protein default `DEFAULTGEP_B = 0` → `penalty_ex
= 0` → all four increments are no-ops. DNA default `DEFAULTGEP_N = 0` →
same. Would have surfaced if a future scoring model or `--exp` override
(distinct from `--ep` which sets matrix offset, not gap extend) routed
non-zero penalty_ex to the group DP. All 265 existing tests still pass.

**Same gap exists in `pairwise_align11`** (`crates/mafft-align/src/profile.rs:910`)
— the flat-penalty 1-vs-1 path used for NW-NS-2 single-merges. It takes
`penalty: f64` only, no `penalty_ex`, so adding the C-equivalent
`mi += fpenalty_ex_i` (with `i < lgth1` guard, mirroring `Galign11.c:1291,1318`)
would require an API change (take `gap: &GapModel` or add a `penalty_ex`
parameter). Currently benign — `pairwise_align11` is reached only by
non-FFT no-constraint 1-vs-1 merges (NW-NS-2 only), and the CLI does
not expose `--exp` (which would set ppenalty_ex). All NW-NS-2 modes
remain byte-identical to C.

### §B.2. Missing FMA `mul_add` outside `match_calc_row`

**Location**: `crates/mafft-align/src/{global,local,genaffine}.rs` —
inner DP loops use plain `a + b` instead of `f64::mul_add`. `profile.rs`
already uses `mul_add` in match_calc_row and in the position-specific
gap candidates (§4 fix).

**C reference**: gcc `-O3 -mfma` (or auto-FMA on `arm64`/`aarch64`) fuses
`a + b * c` into a single-rounding FMA. Two-step Rust arithmetic
produces 1-ULP differences that can flip tie-breaks on flat-landscape
matrices (§4 BL50 root cause).

**Why latent**: All currently-tested matrices have wide enough score
landscapes that the 1-ULP differences don't tip tie-breaks. Will surface
if a future scoring model or substitution-matrix combination produces
flatter score distributions.

**Fix**: Replace `a + b * c` with `b.mul_add(c, a)` wherever it appears
in DP arithmetic. Audit-pass + selective FMA.

### §B.3. `--auto`, `--seed`, `--treein`, `--treeout`, `--memsave`, `--anysymbol`, `--leavegappyregion` — UNIMPLEMENTED

**Location**: `crates/mafft-bin/src/main.rs` — these flags are absent;
the CLI rejects them with "unknown argument". (`--reorder`/`--inputorder`
landed 2026-05-13 — see §AA above.)

**C reference**: `scripts/mafft:237-238` (`--seed`/`--seedtable`),
`scripts/mafft:330-343` (`--anysymbol`), `scripts/mafft:399-403`
(`--treeout`), `scripts/mafft:543-545` (`--memsave`), `scripts/mafft:650`
(`--leavegappyregion`), `scripts/mafft:753-757` (`--treein`),
`scripts/mafft:1290-1340` (`--auto`).

**Severity**: HIGH for feature coverage (users running with these flags
get errors), but does not affect byte-parity of any currently-tested
mode.

**Effort**: Mostly CLI plumbing. `--seed`/`--treein` need engine support
for external guide-tree / seed-alignment input. `--auto` needs a
strategy-selection heuristic mirroring `scripts/mafft:1290-1340`.

### §B.4. `--retree N` for N ≠ 2 byte-untested for INS-i modes

**Location**: `crates/mafft-core/src/engine.rs:410-422` overrides
`retree = 1` for L/G/E/Q/X-INS-i; for FFT-NS-i uses CLI default
`retree = 2`.

**C reference**: `scripts/mafft:86, 142-156`.

**Status**: Logic looks correct (default 2, override to 1 for INS-i),
but no test asserts byte-identity for `--retree 3 --localpair sample` or
similar.

**Effort**: Add a regression test once an output fixture is generated
from C.

### §B.5. HashMap iteration order in `profile_cache`

**Location**: `crates/mafft-core/src/progressive.rs:196, 592, 694`
(`HashMap<Vec<usize>, CachedProfile>`).

**Risk**: Currently safe — the cache is only `.get()`-accessed, never
iterated. If a future refactor adds `.iter()` / `.values()` / `.keys()`
loops over the cache, Rust's randomized iteration order will produce
non-deterministic output.

**Fix**: Add `// DETERMINISM: do not iterate this cache — get-only`
comment, or replace with `BTreeMap` for deterministic iteration order.

### §B.6. Sequential float summation guard in `refinement.rs:1069-1080`

**Location**: `crates/mafft-core/src/refinement.rs:1069-1080` —
intergroup scoring uses explicit sequential summation with a comment
warning against `par_iter().sum()`.

**Status**: Already correct and documented. Listed here as a reminder
during refactoring — removing the guard would cascade into
non-deterministic refinement decisions.

### §B.7. `parttree.rs` (old) `max_by_key` tie-break

**Location**: `crates/mafft-tree/src/parttree.rs:159-161` — uses
`max_by_key` which returns the LAST tied max. C's iteration uses FIRST
tied max.

**Status**: Dead code — the engine uses
`parttree_split::build_parttree_topology` (the resolved §6 path), not
this old `parttree::parttree` function. Reachable only via the deprecated
re-export at `mafft-tree/src/lib.rs:23`. Consider deleting.

---

## §C. Performance / non-correctness items

### §C.1. Per-group gap stripping in progressive alignment

**Status**: Not implemented — we strip columns all-gap across all
sequences globally, not per-group like C's `commongappick()`.

**Effect**: For a 500-column MSA where group1 has 100 residue-containing
columns and group2 has 120, C's per-group strip DP's 100 × 120 = 12,000
cells; ours DP's up to 500 × 500 = 250,000 cells. Functionally correct
(same output) but wastes work.

**Why hard**: per-group stripping breaks when `kept1 != kept2` because
"other" sequences (not in either group) can't follow both cursors
simultaneously. Solving requires porting C's `insertnewgaps()` from
`addfunctions.c` for the strip-restore round-trip.

**Priority**: Low for correctness, medium for performance on large
inputs.

### §C.2. SIMD inner loops

**Status**: SIMD-friendly patterns in `match_score()`,
`pairwise_score()`, `pairwise_identity_distance()` auto-vectorize via
LLVM. The DP fill loops themselves are not SIMD'd (anti-diagonal
restructuring would be required to break the data dependency).

---

## §D. X-INS-i (`--xinsi`) — UNTESTABLE without `contrafold`

Requires Stanford's `CONTRAfold v2.02+` binary, distributed separately
(http://contra.stanford.edu/contrafold/). Not shipped by upstream MAFFT,
not buildable from `mafft-upstream/extensions`. Rust wiring exists
(`engine.rs::XInsi`, `mafft-bin/src/main.rs:--xinsi`) and emits the
correct "contrafold not found" diagnostic when the binary is absent.
End-to-end validation deferred until `contrafold` is installed.

---

## Recommended order of attack

All currently-tested modes are byte-identical to C MAFFT 7.526. The
remaining items are latent / coverage / performance gaps, not active
divergences:

1. **§B.3 missing CLI flags** — `--auto`, `--treein`/`--treeout`,
   `--anysymbol`, `--seed`, `--leavegappyregion`, `--memsave`. Highest
   user-visible impact.
2. **§B.1 `penalty_ex` in `pairwise_align11`** — would need an API
   change (take `&GapModel` or add `penalty_ex` param). Currently
   benign because the CLI doesn't expose `--exp`.
3. **§B.2 FMA `mul_add` in {global,local,genaffine}.rs** — audit pass
   to forestall future 1-ULP tie-break flips on flat-landscape matrices.
4. **§B.4 `--retree N` for N ≠ 2** — add regression test.
5. **§C.1 per-group gap stripping** — performance only, no behavior
   change.
6. **§D X-INS-i (`--xinsi`)** — needs `contrafold` binary to validate.
