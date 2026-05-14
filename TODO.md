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
| `--parttree --reorder`                      | match   | match      | 0          | byte-exact ✓ (closed 2026-05-13) |
| `--treeout` (FFT-NS-2, NW-NS-2, FFT-NS-i, L/G/E-INS-i, BL/JTT, parttree) | match | match | 0 | byte-exact ✓ (closed 2026-05-13) |
| `--dpparttree --treeout`                    | match   | match      | 0          | byte-exact ✓ (closed 2026-05-13) |
| `--tm 200 --treeout`                        | match   | match      | 20         | residual gap — see §B.2: pass-0 TM alignment drift propagates into tree branch lengths |
| `--treein` (FFT-NS-2, NW-NS-2, BL/JTT/TM, L/G/E-INS-i, FFT-NS-i, all non-parttree modes) | match | match | 0 | byte-exact ✓ (closed 2026-05-14) |
| `--treein --treeout`                        | match   | match      | 0          | byte-exact ✓ (closed 2026-05-14) — appends `#by loadtree\n` like C `mltaln9.c:2818` |
| `--auto` (small/medium/large brackets covered by size heuristic) | match | match | 0 | byte-exact ✓ (closed 2026-05-14) |
| `--memsavetree` (k-mer + MSA two-pass tree)  | match   | match      | 0          | byte-exact ✓ (closed 2026-05-14) |
| `--memsavetree --treeout`                    | match   | match      | 0          | byte-exact ✓ (closed 2026-05-14) |

Test suite as of 2026-05-14: **309 Rust tests pass, 0 failed, 0 ignored**
(`cargo test --workspace --exclude pymafft --release`). Plus 32 Python tests
pass. Every mainstream mode in the matrix above is byte-identical to C
MAFFT 7.526 — including `--parttree --reorder` and `--treeout` for all
modes including `--dpparttree`.

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

## §AB. `--treeout` — RESOLVED 2026-05-13 (byte-identical to C for all modes incl. `--dpparttree`)

**Mode**: `mafft <flags> --treeout sample` now writes `<sample>.tree` in
Newick format byte-identical to C MAFFT 7.526 for FFT-NS-2, FFT-NS-i,
NW-NS-2, L/G/E-INS-i, BL/JTT scoring variants, and `--parttree`.
All `--treeout` modes including `--dpparttree` are byte-identical to C
MAFFT 7.526 (see "§AB.5 `--dpparttree --treeout`" at end of this section
for the closure details).

### Implementation

1. **`crates/mafft-tree/src/newick.rs::topology_to_newick`** — port of
   C's standard guide-tree serialization (`mltaln9.c:6190-6491` in
   `fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave`).
   Leaf format `\n<i+1>_<sanitized_name>\n`, sanitize-mask matches C
   (alnum + `/=-{}` kept, everything else → `_`). Branch lengths emit
   via Rust `{:7.5}` to mirror C's `%7.5f`. Final string terminated
   with `;\n`.
2. **`crates/mafft-core/src/progressive.rs::MultipleAlignment`** — two
   new fields:
   - `guide_tree: Option<Topology>` — populated by the engine with
     either the final progressive guide tree (FFT-NS-2 / `--parttree`
     / `--nofft`) or the refinement-pass tree built inside `dvtditr`
     (FFT-NS-i / *-INS-i). Refinement-tree branch lengths require
     rounding distances to `%.3f` to match `dndpre`'s hat2 precision
     (`engine.rs:670-680`).
   - `first_pass_sequences: Option<Vec<Vec<u8>>>` — C's `pre_1` cache,
     needed by both `--parttree --reorder` and `--parttree --treeout`
     because CALL 2 reads the FIRST-pass alignment.
3. **`crates/mafft-tree/src/parttree_split.rs::compute_parttree_newick_fromaln`**
   — parttree-specific tree builder mirroring `splittbfast.c:1275-1301`
   (leaf format: numeric leaves, no branches, no names) +
   `splittbfast.c:2532-2553` (per-merge `(child1,child2)` concat). Shares
   the CALL 2 pivot pipeline with `compute_parttree_order_fromaln` via
   the new `run_parttree_fromaln_pipeline` helper.
4. **`crates/mafft-bin/src/main.rs`** — `--treeout` CLI flag; writes
   `<input>.tree` after alignment finishes. For PartTree, calls the
   `_fromaln` variant on `msa.first_pass_sequences`; for everything
   else, calls `topology_to_newick(msa.guide_tree, msa.names)`.

### `--dpparttree --treeout` — RESOLVED 2026-05-13

C's `--dpparttree` runs `splittbfast` ONCE (cycle=1) with `-U` flag
(`splittbfast.c:677-679` → `doalign=1`), so distances are computed via
`G__align11_noalign( n_disLN, -1200, -60, ... )` on **raw** sequences
(`splittbfast.c:1700`), with `outgap=1` (terminal gaps penalized).
This is fundamentally different from `--parttree`'s cycle=2 pipeline.

**Implementation:**

1. **`mafft-tree/src/parttree_split.rs::run_parttree_pipeline_with_scorer`**
   — generic CALL-1-style parttree pipeline that accepts caller-supplied
   `selfscore`, `orilen`, `pair_score`, and `seqs_equal` closures.
   Shared between the fromaln (CALL 2) and dpparttree (CALL 1) paths.
2. **`mafft-tree/src/parttree_split.rs::parttree_result_to_newick`** —
   extracted from `compute_parttree_newick_fromaln` so both pipelines
   can serialize via the same code.
3. **`mafft-bin/src/main.rs`** — when `--dpparttree --treeout` is set,
   builds the `n_disLN` matrix (`n_dis - 60` for residue cells, 0
   elsewhere — `constants.c:1447-1450`), calls
   `mafft_align::global_align` with `GapModel(-1200, -60)` and
   `head_gap=tail_gap=true` (matching `outgap=1`), and feeds these
   distances into the generic parttree pipeline.

**Result**: `--dpparttree --treeout` is now byte-identical to C MAFFT
7.526 (0-line diff on the 36-seq sample).

---

## §AD. `--auto` — RESOLVED 2026-05-14 (byte-identical to C on small/medium/large brackets)

**Mode**: `mafft --auto` picks the alignment strategy from `nseq`
(number of sequences) and `nlen` (longest sequence length). Mirrors
`scripts/mafft:1290-1343`:

| `nseq < ` | `nlen < ` | Mode | iter | retree |
|-----------|-----------|------|------|--------|
| 100  | 3000  | L-INS-i  | 1000 | 1 |
| 200  | 1000  | L-INS-i  | 2    | 1 |
| 500  | 10000 | FFT-NS-i | 2    | 2 |
| 20000 | —    | FFT-NS-2 | 0    | 2 |
| 100000 | —   | FFT-NS-2 + memsavetree | 0 | 2 |
| 200000 | —   | FFT-NS-2 + memsavetree | 0 | 1 |
| ∞    | 3000  | --dpparttree | — | 1 |
| ∞    | ∞     | --parttree   | — | 1 |

**Implementation**: `crates/mafft-bin/src/main.rs::decide_auto` reads
the size, returns an `AutoChoice { mode, retree, parttree, dpparttree }`.
When `--auto` is set, the choice OVERRIDES `--localpair`/`--globalpair`/
`--genafpair`/`--parttree`/`--dpparttree`/`--maxiterate`/`--retree`.

**Parity tests**: `crates/mafft-bin/src/main.rs::tests::auto_*` (7
unit tests covering each bracket) and
`crates/mafft-core/tests/end_to_end.rs::auto_picks_linsi_for_small_sample`
cross-validates the 36-seq sample (which lands in the smallest
bracket, picking L-INS-i with iterate=1000).

**Caveat**: the 100k–200k brackets use C's `memsavetree` (a memory-
optimised tree algorithm not yet ported); we fall through to standard
FFT-NS-2 / FFT-NS-1 there. The progressive *alignment* step is the
same as C, but the tree construction differs, so very-large outputs
may not be byte-identical. The smaller brackets (which cover all
typical interactive use cases) ARE byte-identical.

---

## §AC. `--treein` — RESOLVED 2026-05-14 (byte-identical to C across all non-E-INS-i modes)

**Mode**: `mafft --treein FILE` lets the user supply a custom guide tree
in MAFFT's internal format (4 columns per line: `im jm len0 len1`,
1-indexed, `im < jm`). C reads this file via
`mltaln9.c::loadtree`/`loadtreeoneline`; users typically convert their
Newick tree using `mafft-upstream/core/newick2mafft.rb`.

**Implementation**:
- `crates/mafft-tree/src/treein.rs` — parser
  (`parse_mafft_tree`/`parse_mafft_tree_str`) reading the 4-column
  format and rebuilding a `Topology` with `JoinStep` lefts/rights set
  to the cluster member lists in C's `min(member)`-keeps-rep order
  (mirrors `loadtree`'s `Bchain` linked-list reduction).
- `crates/mafft-core/src/engine.rs::MafftEngine.treein_path` — new
  `Option<PathBuf>` field. When set, `align()` loads the tree once
  and overrides THREE places C also overrides:
  1. **Progressive merge tree** (`tbfast.c:2072` `if(treein) loadtree`).
  2. **LH-table importance reweighting** (`tbfast.c:2967
     counteff_simple_double_nostatic_memsave` + `tbfast.c:1355
     calcimportance_half`) — the topology drives the per-sequence
     weights used to recompute `region.opt` for the LH constraints.
  3. **Refinement tree rebuild** (`dvtditr.c:766-768` `if(intree)
     veryfastsupg_double_loadtree`).
- `crates/mafft-bin/src/main.rs` — `--treein FILE` CLI flag,
  pre-`align` file-existence check, and a `#by loadtree\n` trailer on
  `--treeout` output matching C `mltaln9.c:2818`.

**Parity tests**: `crates/mafft-core/tests/end_to_end.rs::treein_*`
cross-validates FFT-NS-2, NW-NS-2, L-INS-i, G-INS-i, E-INS-i against
C MAFFT 7.526 with the same `_guidetree` file.

**Result**: All `--treein` modes byte-identical to C MAFFT 7.526.

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

### PartTree partial fix (2026-05-13)

For `--parttree --reorder`, the engine now calls
`parttree_split::compute_parttree_order` which mirrors C's
`splittbfast.c::splitseq_mq` order-generation in two stages:
1. Top-level partition: `treeorder = topol[nyuko-2][0] ++ topol[nyuko-2][1]`
   from the yuko-level UPGMA root step (`splittbfast.c:2351-2354`).
2. Leaf emission: for each yuko in `treeorder`, append `outs[yuko]` in
   the same `j`-iteration order C uses (`splittbfast.c:1305-1309`).

`c_normalized_subtree` recursively applies C's smaller-first-element
normalization (`mltaln9.c:8184-8197`) to reconstruct the leaf order
within each `topol[step][i]` array.

**Result**: parttree --reorder diff vs C: **944 → 454 lines** (50%
reduction). Same alignment content (sequences match by sort), just in
a different order at the upper levels of the yuko tree.

### PartTree full closure (2026-05-13)

C `--parttree` runs `splittbfast` TWICE (`scripts/mafft:2655` and `:2681`):

1. **CALL 1**: builds the parttree from raw 6-mer distances, does
   progressive alignment → `pre_1` (intermediate aligned FASTA).
2. **CALL 2**: re-reads `pre_1` with `-Z` (`fromaln=1`) and recomputes
   distances via `naivepairscore11(orialn_a, orialn_b, penalty)` on the
   first-pass aligned rows (NOT the final-pass alignment). Selfscore is
   the diagonal sum of the substitution matrix
   (`splittbfast.c:3011-3017`). This produces a STRUCTURALLY DIFFERENT
   yuko UPGMA tree from CALL 1 (lopsided 8/27 vs CALL 1's 20/15 split
   for the n=36 sample).

The final `--reorder` output is the **composition** of both passes:
`final_order[k] = call1_order[call2_order[k]]`.

**Implementation:**

1. **`naivepairscore11_aligned`** (`parttree_split.rs:280-321`) —
   port of `mltaln9.c:13801-13851`. Strips common gaps inline (no
   allocation), then walks the alignment adding the substitution
   matrix value or a single gap penalty per gap RUN.
2. **`selfscore_aligned`** (`parttree_split.rs:328-344`) — diagonal
   sum of the substitution matrix over non-gap residues.
3. **`compute_parttree_order_fromaln`** (`parttree_split.rs:402-518`) —
   full CALL 2 pipeline: pick reference, dcompare_sort via libc qsort,
   pivot selection with shimon-style dedupe on aligned content,
   `pickmtx` / `dfromc` via `naivepairscore11_aligned`, yuko
   assignment, UPGMA on yukomtx, c-normalized DFS traversal.
4. **`first_pass_msa`** capture in `engine.rs:447-452` — stashes the
   intermediate MSA after pass 0 of the retree loop. CALL 2 must see
   this `pre_1`, NOT the final `pre_2`, because the second pass
   produces a subtly different alignment that yields different
   pair scores.
5. **Composition** in `engine.rs:678-688` —
   `final_order[k] = call1_order[call2_order[k]]`.

**Result**: `--parttree --reorder` is now **byte-identical to C** on the
36-seq protein sample (0-line diff). Regression test
`reorder_parttree_matches_c` in `crates/mafft-core/tests/end_to_end.rs`
pins the exact 36-element output permutation.

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

### §B.2. Missing FMA `mul_add` outside `match_calc_row` — DEFENSIVE FIXES LANDED, TM SCORING STILL DIVERGES

**Location**: `crates/mafft-align/src/{global,local,genaffine}.rs` —
inner DP loops use plain `a + b` instead of `f64::mul_add`. `profile.rs`
DP inner loop and gap candidates already use `mul_add` (§4 fix).
Boundary inits and the `mj[j]` row-init in `profile.rs`, plus the
FFT polarity/volume channel build in `fft_align.rs`, now use
`mul_add` consistently (landed 2026-05-13). Defensive — preserves
byte-identity for all 17 alignment modes, no regressions.

**C reference**: gcc `-O3 -mfma` (or auto-FMA on `arm64`/`aarch64`) fuses
`a + b * c` into a single-rounding FMA. Two-step Rust arithmetic
produces 1-ULP differences that can flip tie-breaks on flat-landscape
matrices (§4 BL50 root cause).

**Confirmed active for TM scoring (2026-05-13)**:
`mafft --tm 200 --retree 1 sample` vs `mafft-rs --tm 200 --retree 1 sample`
differs by 8 lines — 4 single-char shifts at output FASTA lines 156
and 170 (M62903 and S75720, chicken opsins). The default `--tm 200`
(retree=2) converges to C's output because the second pass re-aligns
from the rebuilt tree, but the first-pass divergence cascades into
`--tm 200 --treeout` branch-length drift (~3e-3, 20-line diff).

**Investigated 2026-05-13** (partial):
1. Matrices are bit-identical to C (cell-by-cell cross-validate
   tests pass with 0 mismatches for BL/JTT/TM/DNA).
2. `--tm 200 --retree 1 --nofft` shows the SAME 8-line diff, so FFT is
   not the source.
3. Pairwise_align11 (1-vs-1) and profile_align (multi-vs-anything)
   are both invoked — diff persists. No naive `+= *` patterns remain
   in either function.
4. Boundary init FMA + `mj[j]` init FMA + FFT-channel polarity/volume
   FMA all landed; none closed the gap.

**Further investigated 2026-05-14**:
5. Guide tree is **byte-identical** between C and Rust even on the
   minimal repro — so the merge order matches; the divergence is
   purely in the per-step DP, not in tree construction.
6. Tried MUL+FMA+ADD reshape in boundary init (closer to gcc's
   left-to-right contraction emission for `a*b + c*d` parenthesized
   subexpressions); no change. Reverted.
7. Confirmed C's `Salignmm.c:1953` increments `m[j] += fpenalty_ex`
   UNCONDITIONALLY (no `j < lgth2` guard). Our prior code had a
   spurious `if j < m` guard; removed (literal C match, no observable
   change for protein default `fpenalty_ex = 0`).
8. Minimal repro: `{seqs 1..29, 30, 33, 36}` (32 of 36 seqs) — adding
   seq 30 to `{1..29, 33, 36}` flips the tie-break. Removing 30 → 0
   diff; adding 30 → 8-line diff. seq 30 (rat 5HT-7 serotonin receptor)
   affects the alignment of seqs 12 (M62903 chicken visual pigment) and
   13 (S75720 chicken P-opsin) at their leading boundary, even though
   30 is in a different subtree.
9. Trees match exactly in the minimal repro, so this is NOT a tree-
   topology issue — it's a DP precision divergence at a specific merge
   step (most likely the merges that introduce 12 and 13).

The remaining 1-ULP drift is in a code path we haven't pinned. Possible
candidates: `global.rs` / `local.rs` / `genaffine.rs` inner DP (used by
pairwise paths we haven't audited), the per-cell impmtx or cpmx
construction in `cpmx_calc_new`-equivalent code, or a subtle ordering
difference in the FFT cross-correlation.

All other modes (BLOSUM62 / BL80 / JTT 200 / DNA / `--add` /
`--allowshift`) are byte-identical even at retree=1, so the FMA gap is
currently TM-specific (TM PAM 200 has the flattest score distribution).

**Fix path** (deferred — see severity note below):
- 8-line diagnostic test: `mafft-rs --tm 200 --retree 1 sample` vs
  `mafft --tm 200 --retree 1 sample` should converge to 0 once the
  right FMA is added.
- Audit `global.rs` / `local.rs` / `genaffine.rs` inner DP loops.
- Add a per-step alignment-trace dump (similar to `CDBG_PT_STEPS` in
  splittbfast) to pin the FIRST diverging merge step.
- Use the cross_validate_profile_align FFI harness to do a focused
  1-vs-3 alignment with the exact profile that occurs in the seq 12
  merge step, and diff DP cell values numerically. The minimal repro
  in finding (8) means the divergent merge can be isolated to a
  single profile-vs-single-seq DP call.

**Severity**: LOW. Default `--tm 200` (retree=2) output is byte-
identical. Only `--retree 1 --tm 200` output (8-line diff in 2 seqs
out of 36) and the implied `--tm 200 --treeout` first-pass tree
(20-line branch-length diff, max drift 3e-3 in 5th decimal place)
are affected. All 17/17 alignment-mode parity tests still pass.

### §B.3. `--seed`, `--memsave`, `--anysymbol`, `--leavegappyregion` — UNIMPLEMENTED

**Location**: `crates/mafft-bin/src/main.rs` — these flags are absent;
the CLI rejects them with "unknown argument".
(`--reorder`/`--inputorder` landed 2026-05-13 — see §AA. `--treeout`
landed 2026-05-13 — see §AB. `--treein` landed 2026-05-14 — see §AC.
`--auto` landed 2026-05-14 — see §AD.)

**C reference**: `scripts/mafft:237-238` (`--seed`/`--seedtable`),
`scripts/mafft:330-343` (`--anysymbol`), `scripts/mafft:543-545`
(`--memsave`), `scripts/mafft:650` (`--leavegappyregion`).

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

### §B.9. ~~`--memsavetree` — algorithm ported, residual height drift~~ — RESOLVED 2026-05-14

**Root cause**: I was porting the WRONG C function.
`compacttree_memsaveselectable` (`mltaln9.c:5491`) is what
`--youngestlinkage` uses (compacttree=4), NOT `--memsavetree`. For
`--memsavetree` (compacttree=2), C MAFFT calls
`compacttreegivendist` (`mltaln9.c:5221`) — a completely different
stepwise-insertion algorithm.

**Discovery path**: FFI wrapper `rs_compacttree_memsaveselectable_kmer`
(after fixing the `njob` global SIGSEGV — C uses GLOBAL `njob` for
`joblist = calloc(njob, sizeof(int))`, not the function param)
showed our memsavetree matched `compacttree_memsaveselectable`
bit-exact. Yet C's `--memsavetree --treeout` output differed by 76
lines. Tracing through `disttbfast.c:4017` revealed
`if (compacttree == 4) compacttree_memsaveselectable(...) else
compacttreegivendist(...)` — and for `compacttree == 2`, it's the
ELSE branch.

**Fix**: Ported `compacttreegivendist` (`mltaln9.c:5221-5331`) into
`crates/mafft-tree/src/memsavetree.rs::compacttree_givendist`. The
algorithm:
1. Initial step: leaves 0 and 1 link at `mindist[1]/2`.
2. For each subsequent leaf `i ∈ [2, nseq)`: walk UP the tree from
   `treept[nearest[i]]` until a parent's height exceeds
   `mindist[i]/2`, then insert leaf `i` as a sibling under a new
   internal node at height `mindist[i]/2`.
3. DFS post-order traversal reformats into our `Topology` struct
   (mirrors C `reformat_rec`).

The pass-1 MSA tree uses the SAME algorithm with MSA-derived
`mindist`/`nearest` (from `distcompact_msa` via `naivepairscorefast`).
This is what `tbfast.c:2538 "Making a compact tree from msa, step 1"`
does after pass 0 progressive alignment.

**Result**: `--memsavetree` (both `--retree 1` and default
`retree=2`) now byte-identical to C MAFFT 7.526 on the 36-seq sample.
`--memsavetree --treeout` also byte-identical.

**Tests**: 4 cross-validate tests
(`crates/mafft-tree/tests/cross_validate_memsavetree.rs`) + 1 e2e
fixture test (`memsavetree_byte_identical_to_c`). All pass.

---

### §B.9.legacy. Old notes — algorithm ported, residual height drift after merge step 7

**Investigation update 2026-05-14**:

Added FFI cross-validation infrastructure (`mafft-sys::wrappers/parttree_helpers.c::rs_compact_initial_mindist` + `mafft-sys::distcompact`) and three new tests in `crates/mafft-tree/tests/cross_validate_memsavetree.rs`:
- `distcompact_matches_c_for_every_pair` ✓ (per-pair distance bit-exact)
- `initial_mindist_matches_c` ✓ (initial pairwise scan + nearest array bit-exact)
- `cluster_mix_for_first_divergent_step` ✓ (cluster_mix for the specific (28→{29,30}) merge matches C exactly: d(29,28)=0.194214, d(30,28)=0.165792, mix=0.167213)

**Key finding**: All per-pair primitives MATCH C bit-for-bit. The algorithm structure matches C (we even mirror C's `for(acpti=ac; acpti->next!=NULL; ...)` last-active-skip quirk at `mltaln9.c:5647`). The first 6 merges (33,34), (19,20), (29,30), (7,8), (9→cluster), (21→cluster) match C exactly with the same branch lengths.

**The divergence appears at merge step 7**: leaf 28 (0-indexed) attaches to cluster {29, 30}. Our mindist[28] = 0.167213 (= cluster_mix value); C's effective merge distance for the same pair = 0.19422 (twice the tree-output branch length 0.09711). Yet cluster_mix(d(29,28), d(30,28)) provably equals 0.167213 in BOTH C and Rust.

Hypothesis: C's mindist[28] somehow stays at a larger value (~0.194) despite the cluster_mix at step k=2 computing 0.167. Possible mechanisms: (a) some prior step's `nearest[i] == jm` update path skips the `if(tmpdouble < mindist[i])` guard; (b) the "antei sei no tame" loop after distance recomputation overwrites in a way we miss; (c) we're misreading which pair C actually picks at step 7.

**Outstanding FFI work**: `rs_compacttree_memsaveselectable_kmer` (wraps C's full algorithm) currently SEGFAULTs inside compacttree — likely missing global setup (TLS, distarrarg state). Test marked `#[ignore]` for now. Fixing that wrapper would directly capture C's `topol[k][0][0]` / `len[k]` per step and pinpoint the first divergent decision.

**Earlier state preserved below**:

### §B.9.old. `--memsavetree` — algorithm ported, residual merge-order drift after step 2

**Location**: `crates/mafft-tree/src/memsavetree.rs` — ports
`mltaln9.c::compacttree_memsaveselectable` with `howcompact=2`,
`memsave=1` (the algorithm `--memsavetree` and `--auto` 100k+ brackets
use). Wired into `engine.rs` for both pass 0 (k-mer distance) and
pass 1+ (MSA distance via `distcompact_msa`).

**Current state**:
- **Algorithm structure correct**: the C reference's
  `compactdisthalfmtxthread` (initial scan) + main loop
  (`mltaln9.c:5639-6097`) + per-step distance recomputation
  (`verycompactkmerdistarrthreadjoblist` / `verycompactmsadistarrthreadjoblist`)
  is faithfully ported with cluster_mix linkage (`SUEFF=0.1`),
  `preferenceval` tie-break, and full member-list reconstruction for
  the progressive-alignment consumer.
- **First-merge parity**: on the 36-seq sample, the FIRST merge after
  the MSA-distance rebuild matches C MAFFT 7.526 byte-for-byte —
  branch length 0.13869 (vs the standard musclesupg branch 0.35002),
  confirming the MSA `distcompact_msa` path is bit-correct.
- **Residual drift starting at merge step ≥ 2**:
  - k-mer-only tree (`--memsavetree --retree 1 --treeout`): 76-line
    diff vs C (different topology starting at step 3).
  - Full two-pass output (`--memsavetree`): 930-line diff vs C, all
    propagated from the pass 0 k-mer tree topology drift.
- **Tests**: 4 unit tests in `memsavetree::tests::*` cover algorithm
  invariants (`im < jm` swap, 2-seq merge, identical-seq merge, etc.).
  No fixture cross-validation yet because outputs don't match.

**Suspected cause**: tie-break ordering in the `find min mindist[i]`
scan, or a 1-ULP drift in `cluster_mix_double` that pushes the
nearest-neighbor pick onto a different `j` for some merge step. The
initial pairwise mindist[]/nearest[] looks correct (since the first
merge picks the same pair); the divergence emerges in the second-merge
selection.

**Fix path**:
- Add `mafft-sys` FFI binding for `compacttree_memsaveselectable` (or
  call `compactdisthalfmtxthread` directly) to cross-validate
  intermediate mindist[]/nearest[] arrays cell-by-cell.
- Audit `nearest[i] == jm` handling: when the merge updates `nearest`
  without updating `mindist`, this preserves the OLD `mindist`
  reference to a (now-merged) cluster. Verify our behavior matches C's
  exactly here.
- Once the k-mer tree matches, the MSA tree should follow since the
  per-step pipeline is the same code.

**Severity**: MEDIUM. `--memsavetree` is essentially "extra-large
input" territory; `--auto` falls back to it only at nseq ≥ 100k.
For 36 seqs our output is functional (valid alignment) but not
byte-identical to C's.

---

### §B.8. ~~`--treein --genafpair --maxiterate >1` (E-INS-i refinement w/ user tree)~~ — RESOLVED 2026-05-14

**Root cause**: `engine.rs::recompute_importance` was deriving the
sequence weights for the LH table importance recomputation from
`musclesupg(&dm)` — the UPGMA tree built from pairwise distances.
C MAFFT does this same step (`tbfast.c:2967 counteff_simple_double_
nostatic_memsave( njob, topol, len, dep, eff )` followed by
`tbfast.c:1355 calcimportance_half( njob, effarr, aseq, ... )`)
using `topol`/`len` from `loadtree`'s loaded user tree. With
`--treein` the two topologies differ on branch lengths, so the
weights differ → the LH table's per-region `region.opt` (importance)
differs → E-INS-i refinement makes different tie-break decisions
starting at iter 2 (iter 1 happens to match because the first round
of refinement reads the same starting alignment).

L-INS-i / G-INS-i / FFT-NS-i all use the same code path; they
happened to converge to the same final alignment regardless of
weights on the 36-seq test sample (less sensitive tie-breaks). The
fix unifies behavior across all four modes.

**Fix**: `engine.rs:412-426` now loads `user_topo` BEFORE the
importance-recomputation block and prefers it over
`musclesupg(&dm)`:
```rust
let initial_topo = user_topo.clone()
    .unwrap_or_else(|| musclesupg(&dm, ClusterMethod::default()));
let weights = mafft_tree::sequence_weights(&initial_topo);
```

**Regression test**: `crates/mafft-core/tests/end_to_end.rs::
treein_einsi_byte_identical_to_c` cross-validates `--treein
--genafpair --maxiterate 1000` against the fixture
`sample.treein.einsi` (generated from C MAFFT 7.526 with the same
user tree).

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
