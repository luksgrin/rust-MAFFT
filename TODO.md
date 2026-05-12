# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION,
> not just at convergence. A run that matches C's final width by coincidence
> while producing different intermediate alignments is NOT parity — it means
> our accept trajectory differs from C's, and the match is fragile
> (input-dependent).

## Current parity matrix (36-seq protein sample, mafft-upstream/test/sample)

Verified 2026-05-11 by running `target/release/mafft-rs <args> sample` against
`mafft-upstream/scripts/mafft <args> sample` on the upstream submodule
(MAFFT 7.526) and diffing the FASTA outputs.

| Mode / flags                                | C width | Rust width | diff lines | Status |
|---------------------------------------------|---------|------------|------------|--------|
| FFT-NS-2 (default)                          | 717     | 717        | 0          | byte-exact ✓ |
| NW-NS-2 (`--nofft`)                         | 717     | 717        | 0          | byte-exact ✓ |
| FFT-NS-i (`--maxiterate 100`)               | 721     | 721        | 0          | byte-exact ✓ |
| L-INS-1 (`--localpair --maxiterate 0`)      | 719     | 719        | 0          | byte-exact ✓ |
| G-INS-1 (`--globalpair --maxiterate 0`)     | 746     | 746        | 0          | byte-exact ✓ |
| E-INS-1 (`--genafpair --maxiterate 0`)      | 729     | 729        | 0          | byte-exact ✓ |
| L-INS-i (`linsi` / `--localpair --maxiterate 1000`) | 731 | 731 | 0       | byte-exact ✓ |
| G-INS-i (`ginsi` / `--globalpair --maxiterate 1000`)| 746 | 746 | 0       | byte-exact ✓ |
| E-INS-i (`einsi` / `--genafpair --maxiterate 1000`) | 729 | 729 | 0       | byte-exact ✓ |
| BL30/45/50/62/80 FFT (`--bl N`)             | match   | match      | 0          | byte-exact ✓ |
| BL80 NW (`--bl 80 --nofft`)                 | 712     | 712        | 0          | byte-exact ✓ |
| JTT 100/200 FFT (`--jtt N`)                 | match   | match      | 0          | byte-exact ✓ |
| TM 100/200 NW + FFT (`--tm N`)              | match   | match      | 0          | byte-exact ✓ |
| PartTree (`--parttree`, `--dpparttree`)     | 752     | 752        | 0          | byte-exact ✓ |
| `--add` / `--add --nofft` / `--add --keeplength` | match | match    | 0          | byte-exact ✓ |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 62 (case)  | byte-exact mod case ✓ |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 62 (case)  | byte-exact mod case ✓ (needs `mxscarnamod`) |
| `--allowshift --globalpair sample`          | 1060    | 1018       | many       | partial — warp DP ported; 42-col gap traced to per-pair re-align bug (see §A) |

Test suite as of 2026-05-12: **269 Rust tests pass, 0 ignored** (`cargo test
--workspace --exclude pymafft --release`). Plus 32 Python tests pass.

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

## §A. `--allowshift` — PARTIAL (42-col gap remains)

**Mode**: `mafft --allowshift --globalpair --maxiterate 0 sample` produces
Rust width **1018** vs C width **1060**. All other byte-identity modes still
pass.

### Mechanism (C MAFFT 7.526)

`--allowshift` triggers four cooperating effects:

1. **`unalignlevel = 0.8`** (C's `specificityconsideration`).
2. **Per-step dynamic matrix scaling** in progressive merge
   (`disttbfast.c:2304`, `tbfast.c:1440`, `mltaln9.c::makedynamicmtx`):
   each merge step scales the substitution matrix by `min(0,
   distfromtip - unalignlevel) * 600`.
3. **Per-pair dynamic re-alignment** in `pairlocalalign.c:2197-2228`:
   after initial pairwise alignment, if `0.5 * dist - unalign_level < 0`,
   re-run the pairwise DP with a dynamic matrix.
4. **Warp DP** activated by `penalty_shift_factor = spfactor = 2.0` (<10
   triggers `trywarp = 1` in `constants.c:277-278`). The recurrence in
   `Galign11.c:780-870` and `Salignmm.c:1957-2020` allows cell (i,j) to
   "warp" back to an anchor (i',j') with cost `penalty_shift +
   penalty_ex * Manhattan(i-i', j-j')`.
5. **Script param zeroing** (`scripts/mafft:1469-1473`): when
   `unalignlevel > 0`, `laof = lexp = pgaof = pgexp = 0`, making the
   pairwise pscores larger (no per-cell offset/extend penalty).

### What's done

| Effect | Rust impl | Status |
|--------|-----------|--------|
| `unalign_level` field on engine | `engine.rs:65,151` | ✓ |
| CLI `--allowshift`, `--unalignlevel` | `main.rs:48-58` | ✓ |
| Per-step matrix scaling | `progressive.rs::progressive_align_full`, `dist2offset`, `make_dynamic_matrix` | ✓ |
| Per-pair dynamic re-alignment | `constraints.rs::build_homology_table_with_unalign` | ✓ |
| Warp DP in pair phase (`global_align`) | `global.rs:120-265` | ✓ byte-identical for single pairs |
| Warp DP in group phase (`profile_align_imp_with_boundary`) | `profile.rs:587-737` | ✓ byte-identical for single 1-vs-1 |
| Script param zeroing (lexp=laof=pgexp=pgaof=0) | NOT applied | ✗ blows width to 1412+ when applied |

### Verification

- `rust_global_align_matches_c_g__align11_warp` (cross_validate_constrained_align.rs)
  — single-pair warp DP byte-identical to C's `G__align11` with
  `penalty_shift_factor = 2.0`.
- `profile_align_imp_warp_matches_c_a__align` (cross_validate_profile_align.rs)
  — profile-DP warp byte-identical to C's `A__align` for a single 1-vs-1
  group merge.

### Investigation 2026-05-12: warp DP cleared, bug isolated to profile_align_imp + impmtx + shifted matrix

Bisecting the 36-seq input by size shows the divergence first appears at
n=11 (no diff for n≤10). At n=14, only ONE pair out of 91 diverges:
(U22180, M62903) — rat opsin vs chicken visual pigment.

What was initially diagnosed as a "global_align warp DP bug under shifted
matrices" turned out to be a **test setup bug**: the C-side test wrote
`penalty=-1199` directly but then called `constants()` AGAIN (to refresh
`trywarp` from `penalty_shift_factor`), and `constants()` re-derived
penalty from `ppenalty` (which was `NOTSPECIFIED` → default `-1530` →
`penalty=-917`). With `ppenalty=-2000` explicitly written before the
second `constants()` call (mirroring `--op 2.00`), the C path uses the
correct penalty=-1199 and Rust and C produce **byte-identical** output
for the re-aligned pair (score 44908.394 in both).

**Tests added (all pass, 2026-05-12)**:
- `rust_global_align_matches_c_g__align11_warp_k03494_m92036` — pair
  (K03494, M92036) with warp, byte-identical.
- `rust_global_align_matches_c_g__align11_warp_u22180_m62903` — pair
  (U22180, M62903) initial, byte-identical.
- `rust_global_align_realign_matches_c_g__align11_warp_u22180_m62903` —
  pair (U22180, M62903) RE-ALIGN with delta=-161, byte-identical (was
  `#[ignore]` thought to surface a bug, now passing after test fix).
- `profile_align_imp_warp_shifted_matrix_matches_c_a__align` — profile
  DP with shifted matrix, byte-identical.

### Bug is NOT in warp DP; it's downstream

`global_align` and `profile_align_imp` warp DPs are both verified
byte-identical to C for the shifted-matrix case. The remaining 42-col
pipeline gap on n=36 must come from one of:

1. **`profile_align_imp` + impmtx interaction under a shifted matrix**.
   The pipeline calls `profile_align_imp(prof1, prof2, dyn_matrix, gap,
   true, true, Some(impmtx))`. None of the unit tests cover this exact
   combination yet.
2. **`build_imp_matrix` divergence vs C's `imp_match_init_strict`** for
   the post-warp pair alignment.
3. **Per-step `makedynamicmtx` delta** at the GROUP DP level (vs per-pair
   delta for the constraint phase). Both should be the same formula but
   they're computed at different points.

### Engine-side fix (NOT applied)

`scripts/mafft:1469-1473` zeros `lexp = laof = pgexp = pgaof = 0` when
`unalignlevel > 0`. Applying that in `engine.rs:304-318` closes the
pair-phase param mismatch BUT amplifies the (1)/(2)/(3) downstream bug
above — width gap balloons from 42 → 411 cols on n=36. Until the
downstream bug is fixed, leaving the zeroing OFF is the better trade.

### Other observations

- The `§B.1` `mi += penalty_ex` fix incidentally closed n≤13 divergences.
- The 42-col gap on n=36 cascades from the impmtx-side bug, not the
  warp DP itself.

### Priority: Low

`--allowshift` is rarely used. Current implementation makes the flag
take meaningful effect (746 → 987, vs C's 1029) without regressing any
mainstream mode.

### Trajectory (for context)

C output for n=36 with --allowshift: 1060 cols.

Rust progression on n=36:
- 746 (no effect, only CLI flag wired)
- 809 (per-step matrix scaling, int matrices)
- 957 (+ per-pair re-align, int matrices)
- 993 (after f64 DP migration)
- 987 (after warp DP port — single-pair byte-identical for some pairs)
- 1018 (after `§B.1` `mi += penalty_ex` fix — closed n=6/10/13 divergences)
- **Target: 1060** (42 cols away; root cause is the re-align bug in §A.2)

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

### §B.3. `--auto`, `--seed`, `--treein`, `--treeout`, `--memsave`, `--reorder`, `--inputorder`, `--anysymbol`, `--leavegappyregion` — UNIMPLEMENTED

**Location**: `crates/mafft-bin/src/main.rs` — these flags are absent;
the CLI rejects them with "unknown argument".

**C reference**: `scripts/mafft:237-238` (`--seed`/`--seedtable`),
`scripts/mafft:330-343` (`--reorder`, `--inputorder`, `--anysymbol`),
`scripts/mafft:399-403` (`--treeout`), `scripts/mafft:543-545`
(`--memsave`), `scripts/mafft:650` (`--leavegappyregion`),
`scripts/mafft:753-757` (`--treein`), `scripts/mafft:1290-1340`
(`--auto`).

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

1. **§A `--allowshift` 42-col gap** — needs multi-pair byte-identity
   test of constraint table + per-step trace dump comparison vs C to
   pin the divergence source. Low priority unless someone actually uses
   the flag.
2. **§B.1 `penalty_ex` accumulation** — add a regression test for
   `--ep 0.5` on the 36-seq sample, then add the missing `mi += f_ext`
   / `mj[j] += f_ext` lines.
3. **§B.3 missing CLI flags** — `--auto` and `--treein` / `--treeout`
   are the highest-value additions.
4. **§C.1 per-group gap stripping** — performance only, no behavior
   change.
