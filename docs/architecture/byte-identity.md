# Byte-identity

rust-MAFFT's output matches C MAFFT 7.526 byte-for-byte across BAliBASE
3 (1930/1930 fixtures) for every supported mode, *when both are built on
the same platform*. This page explains why that's the bar, what it cost
to hit it, and why the "same platform" qualifier is not a hedge but a
property of the C reference itself.

## Why hold to byte-identity

MAFFT is the *de facto* MSA reference. Every downstream tool and
pipeline — phylogenetics, structure prediction, comparative genomics —
treats its output as ground truth. A Rust port that's "close" forces
every downstream pipeline to choose:

1. **Re-validate the rewrite** on their own data, accepting that
   alignment differences may cascade through downstream analyses.
2. **Maintain a parallel C dependency** so they can verify rust-MAFFT
   isn't introducing drift.

Both are non-starters for adoption. Byte-identity removes the choice:
the Rust binary is a drop-in replacement, no warm-up trust period
required.

## Parity status

### Modes covered

**Production-ready.** `mafft-rs` produces **byte-identical output to C MAFFT 7.526** across:

- All progressive modes: FFT-NS-2, NW-NS-2, FFT-NS-i, L-INS-1/i, G-INS-1/i, E-INS-1/i
- All scoring matrices: BLOSUM 30/45/50/62/80, JTT 100/200, TM 100/200 (FFT + NW + `--retree 1`)
- Tree variants: `--parttree`, `--dpparttree`, `--memsavetree`, `--youngestlinkage`, `--averagelinkage`, `--minimumlinkage`, `--mixedlinkage`, `--pileup`
- Iteration strategies: `--bestfirst`, `--simplehillclimbing`, `--skipiterate F` (both large- and small-F branches), `--oneiteration`
- Add-sequence modes: `--add`, `--add --reorder`, `--add --maxiterate`, `--keeplength`, `--mapout` / `--compactmapout` (the map and the added-sequence names), `--addfragments` and `--addfragments --keeplength` (`addsingle` on the script's `distance=multi` branch)
- Auxiliary outputs: `--distout`, `--scoreout`, `--nodeout`, `--mapout` / `--compactmapout`, `--treeout`
- DNA strand auto-detection: `--adjustdirection` (k-mer) and `--adjustdirectionaccurately` (DP)
- Refinement controls: `--allowshift`, `--unalignlevel`, `--maxambiguous`, `--minimumweight`, `--nwildcard`, `--leavegappyregion` / `--legacygappenalty`
- Memory strategy: `--memsave` / `--nomemsave` (Hirschberg DP), `--memsavetree`
- Seed alignment: `--seed FILE`, `--seedtable FILE` (single + multi-file)
- Fine-grained gap penalties: `--op`, `--ep`, `--exp`, `--shiftpenalty`, `--lop`, `--lep`, `--lexp`, `--gop`, `--gep`, `--gexp`
- Symlink dispatch via `argv[0]`: `linsi`, `ginsi`, `einsi`, `fftns`, `fftnsi`, `nwns`, `nwnsi`, `qinsi`, `xinsi`

Every progressive merge step matches in score and width; every refinement iteration converges to C's exact alignment.

### Parity matrix (36-seq protein sample, `mafft-upstream/test/sample`)

| Mode / flags                                | C width | Rust width | diff | Status |
|---------------------------------------------|---------|------------|------|--------|
| FFT-NS-2 (default)                          | 717     | 717        | 0    | ✓ byte-exact |
| NW-NS-2 (`--nofft`)                         | 717     | 717        | 0    | ✓ byte-exact |
| FFT-NS-i (`--maxiterate 100`)               | 721     | 721        | 0    | ✓ byte-exact |
| L-INS-1 (`--localpair --maxiterate 0`)      | 719     | 719        | 0    | ✓ byte-exact |
| L-INS-i (`linsi`)                           | 735     | 735        | 0    | ✓ byte-exact |
| G-INS-1 (`--globalpair --maxiterate 0`)     | 746     | 746        | 0    | ✓ byte-exact |
| G-INS-i (`ginsi`)                           | 737     | 737        | 0    | ✓ byte-exact |
| E-INS-1 (`--genafpair --maxiterate 0`)      | 729     | 729        | 0    | ✓ byte-exact |
| E-INS-i (`einsi`)                           | 729     | 729        | 0    | ✓ byte-exact |
| BLOSUM 30 / 45 / 50 / 62 / 80 (FFT + NW)    | match   | match      | 0    | ✓ byte-exact |
| JTT 100 / 200, TM 100 / 200 (FFT + NW)      | match   | match      | 0    | ✓ byte-exact |
| `--parttree`, `--dpparttree`, `--parttree --nofft`, `--parttree --reorder` | 752 | 752 | 0 | ✓ byte-exact |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6) | match | match | 0 | ✓ byte-exact |
| `--add --reorder` (30+6 and the #40 reproducer) | match | match | 0 | ✓ byte-exact. `f2cl -r` applies `topolorderz` of disttbfast's second guide tree (`-E 2`), not the k-tuple tree that drove the merge. On these inputs that order is input order with sequences 11 and 12 swapped, and the second cycle does not move residues |
| `--add --maxiterate 2` (30+6 and the #40 reproducer) | match | match | 0 | ✓ byte-exact. `dndpre -M 2` then `dvtditr -K nadd -B` (no segmentation; a side is skipped when it overlaps the existing sequences but is not exactly that set). No `--thread` stays on `-C 0` |
| `--add --compactmapout` (30+6 and the #40 reproducer) | match | match | 0 | ✓ byte-exact, including added-sequence names. The flag implies `--keeplength`; names that lost residues are prefixed with `Nins:` and the first insertion run |
| `--addfragments` / `--addfragments --keeplength` (30+6 and the #40 reproducer) | 871 / 595 | 871 / 595 | 0 | ✓ byte-exact. Default `--addfragments` is `distance=multi` (`scripts/mafft:2628`): `dndpre -M 2` on the existing rows, `pairlocalalign -Y` (`L__align11`) for existing-to-fragment distances and `hat3`, then `addsingle -d -O -l 2.7` (`importance = opt * 0.5`, one constrained `A__align`, new gaps `=`). `--keeplength` adds `-Y`. The large-input 6-mer fallback is not this path. `--addfragments --localpair` is not pinned: C reads off the end of a buffer on some of those inputs |
| `--add` adversarial (16-existing + 3-added with random insertions) | match | match | 0 | ✓ byte-exact |
| `--pileup`, `--youngestlinkage`, `--averagelinkage`, `--minimumlinkage` | match | match | 0 | ✓ byte-exact |
| `--bestfirst` (FFT-NS-i / L-INS-i / G-INS-i)  | match | match      | 0    | ✓ byte-exact |
| `--skipiterate F` (F ∈ {0.05…1.0})         | match   | match      | 0    | ✓ byte-exact |
| `--oneiteration` (FFT-NS-2 + FFT-NS-i)      | match   | match      | 0    | ✓ byte-exact |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0    | ✓ byte-exact |
| `--allowshift --globalpair --maxiterate 1000` | 1060  | 1060       | 0    | ✓ byte-exact (all 386 BBaliBase 3) |
| `--adjustdirection`, `--adjustdirectionaccurately` (mixed forward/RC DNA) | match | match | 0 | ✓ byte-exact |
| `--memsavetree`, `--memsavetree --treeout`  | match   | match      | 0    | ✓ byte-exact |
| `--auto` (small/medium/large brackets)      | match   | match      | 0    | ✓ byte-exact |
| `--seed FILE`, `--seedtable FILE`           | match   | match      | 0    | ✓ byte-exact |
| `--treein`, `--treeout`, `--treein --treeout` | match | match      | 0    | ✓ byte-exact |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0    | ✓ byte-exact |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 0‡   | ✓ byte-exact (needs `mxscarnamod`) |
| `--nuc`, `--amino` (forced sequence type)   | match   | match      | 0    | ✓ byte-exact |

‡ Q-INS-i requires the external `mxscarnamod` binary.

### DNA and input-handling parity (issue #1 integration)

| Input / flags | Status |
|---------------|--------|
| 120 × 1.4 kb CDS (`mtb_cds_120x1400.fa`), FFT-NS-i `--retree 2 --maxiterate 2` (1742 columns) | ✓ byte-exact |
| Same input, `--thread 1 --maxiterate 1/2/3/1000` (C's `athread` refinement rules) | ✓ byte-exact |
| First 8 sequences of the above, `--retree 2 --maxiterate 1000` (trailing `searchAnchors` region) | ✓ byte-exact |
| panaroo gene clusters (5 × 4 seqs, ~300–2300 bp), `--auto --adjustdirection --thread 1 --nuc` | ✓ byte-exact |
| 3 × 3.6–5.0 kb CDS (`dna_01715.fa`, issue #43), FFT-NS-1 / FFT-NS-2 / FFT-NS-i / `--auto --adjustdirection --thread 1 --nuc` | ✓ byte-exact |
| DNA L-INS-1 / `--auto` on synthetic clusters (pair-phase gap scale ×3, `dndpre` offset 220) | ✓ 60/60 and 30/30 |
| Two-sequence input under `--maxiterate 1000` (C refines a pair exactly once, unweighted) | ✓ byte-exact |
| Input read path corpus (22 files × `default` / `--anysymbol` / `--localpair --maxiterate 0` / `--nuc` / `--add` modes, 71 cells incl. C's exit-1 cases) | ✓ 70/71 (see [What's still out of scope](index.md#whats-still-out-of-scope)) |

`--thread 1` and no `--thread` are *different* C code paths (`tditeration.c:1433`
selects `athread` on `nthread > 0`), with different branch order and
convergence rules; `mafft-rs` reproduces both byte-for-byte. C's own
`--thread N` for N ≥ 2 is nondeterministic, so byte-identity there is not
defined; `mafft-rs` stays deterministic at every thread count.

Residue case matches C MAFFT: nucleotide output is lowercased and protein
output uppercased at read time (C `io.c:1462-1467`, `io.c:1755`), driven by
the detected type or by `--nuc` / `--amino`; `--anysymbol` / `--preservecase`
keep the input's own case.

### PartTree on 50 or more sequences (issue #38)

Below 50 sequences `splittbfast` keeps every pivot (`tokyoripara = 0`
when `picksize > njob`), so the 36-sequence sample exercises only one
level of `splitseq_mq`. From 50 sequences on, redundant pivots are
filtered, yuko groups hold distinct sequences, `splitseq_mq` recurses
into them, and for `nin > 50` the pivot draw consumes libc `rand()`. The
recursive port (`mafft-tree/src/parttree_splitseq.rs`) reproduces all of
it, including the macOS `rand()` sequence and the macOS `qsort` order of
tied distances. Checked against the in-tree C build on the scaling corpus
(`benchmarks/scaling_corpus.py`):

| Input | Flags | Status |
|-------|-------|--------|
| 50 Response_reg proteins cut to 10 % (`tests/fixtures/parttree_respreg50.fa`) | `--parttree`, `--parttree --reorder`, `--parttree --nofft`, `--dpparttree`, `--treeout` | ✓ byte-exact |
| 10 Pfam families × n = 100 (full-length UniProt) | `--parttree`, `--parttree --reorder`, `--parttree --nofft`, `--dpparttree` | ✓ 40/40 |
| 10 Pfam families × n = 250 | same four | ✓ 38/40 |
| 16S rRNA, n = 100 / 250 | `--parttree --nuc` (± `--reorder`, `--nofft` at n = 100) | ✓ 4/4 |

The two n = 250 residuals are `--dpparttree` on PF00069 and PF00400.
There the merged profile grows past 10000 columns and C's `pairalign`
switches to `MSalignmm` for the remaining merges
(`splittbfast.c:1020-1028`; the last 3 and 13 merges respectively);
`mafft-rs` keeps `A__align`, so those final merges differ (PF00400: same
width, different gap placement; PF00069: 4 columns narrower). Every
earlier merge matches in width.

### BBaliBase 3 parity sweep

Sweep across BBaliBase 3 (RV11–RV50, 386 protein test sets) against MAFFT 7.526:

| Mode | Flags | Match | Total % |
|------|-------|-------|---------|
| FFT-NS-i | `--maxiterate 1000` | **386/386** | **100.0 %** |
| L-INS-i  | `--localpair --maxiterate 1000` | **386/386** | **100.0 %** |
| G-INS-i  | `--globalpair --maxiterate 1000` | **386/386** | **100.0 %** |
| E-INS-i  | `--genafpair --maxiterate 1000` | **386/386** | **100.0 %** |
| G-INS-i + allowshift | `--globalpair --maxiterate 1000 --allowshift` | **386/386** | **100.0 %** |

Total: **1930/1930 alignments byte-identical** to C MAFFT 7.526 across the BBaliBase 3 protein corpus.

## Byte-identity versus alignment scores

The usual way to judge a multiple aligner is to score its output against
curated reference alignments: the sum-of-pairs score (SP, the fraction of
reference residue pairs the test alignment also aligns) and the
total-column score (TC, the fraction of reference columns reproduced
exactly), usually restricted to the reference's structurally reliable
core blocks, as in BAliBASE. Consistency measures such as T-Coffee's TCS
play a similar role without a reference. These are the right instruments
for comparing *alignment methods*: they tell you whether one method
recovers structure better than another.

They are the wrong instrument for verifying a *port*. A port has a single
correct answer, the output of the program it reimplements, and the
question is whether it reproduces that answer, not whether its answer is
about as good. The two criteria are related in one direction only:

- **Byte-identity ⇒ identical scores.** SP and TC are functions of the
  alignment; equal bytes give equal SP, TC, TCS, or any other score.
- **Identical (or near-identical) scores ⇏ identical alignments.** SP and
  TC only look at core-block columns and summarise them as a mean, so two
  alignments that place residues differently, or even have different
  widths, can score the same.

Two concrete cases from this project show the gap.

**One rounding, 26 columns.** C MAFFT built in-tree on arm64 (clang fuses
`a*b + c` into a single-rounding FMA) and on x86-64 (gcc does not) give
different alignments of the 36-sequence `mafft-upstream/test/sample` under
`--bl 50 --retree 2 --maxiterate 0`: 712 and 738 columns. Both are the
same C source, and `mafft-rs` reproduces each on its own platform
(`--features fp-contract-none` gives the x86-64 output on any host). The
sample has no reference alignment, so the two are compared with each
other and with MAFFT's own objective:

| | arm64 C (712 cols) vs x86-64 C (738 cols) |
|---|---|
| `diff` lines (`<` + `>`) | 836 |
| Columns identical in both | 636 (86–89 % of either) |
| Aligned residue pairs shared | 196 992 of 206 884 / 204 781 (95–96 %) |
| MAFFT `--scoreout` objective (BLOSUM50) | 68 058.4 vs 68 410.6 (0.5 %) |
| MAFFT `score` utility (BLOSUM62) | 677.64 vs 680.78 (0.5 %) |

On BAliBASE the same pair of C builds can be compared against the
references. In FFT-NS-2, 8 of the 386 sets give different bytes under
the two contraction policies; on all eight SP is identical to three
decimals, and TC is identical on seven and moves by 0.01 on the eighth:

| Set | Width (arm64 / x86-64) | `diff` lines (`<` + `>`) | SP (arm64 / x86-64) | TC (arm64 / x86-64) |
|---|---|---:|---|---|
| BBS11018 | 869 / 869 | 2 | 0.704 / 0.704 | 0.410 / 0.410 |
| BB20002 | 2435 / 2435 | 12 | 0.441 / 0.441 | 0.000 / 0.000 |
| BB20018 | 476 / 476 | 102 | 0.965 / 0.965 | 0.760 / 0.770 |
| BB20039 | 663 / 663 | 2 | 0.944 / 0.944 | 0.340 / 0.340 |
| BB30028 | 643 / 643 | 2 | 0.888 / 0.888 | 0.320 / 0.320 |
| BB40036 | 846 / 846 | 4 | 0.968 / 0.968 | 0.690 / 0.690 |
| BB40048 | 1995 / 1995 | 2 | 0.917 / 0.917 | 0.750 / 0.750 |
| BB40049 | 5681 / 5681 | 8 | 0.957 / 0.957 | 0.460 / 0.460 |

A score-based acceptance test would call every one of these "equivalent".
A port that drifted like this would differ from the C it claims to
reproduce on 2 % of BAliBASE without SP noticing at all.

**Tied traces.** The two residuals in the README's
[Known limitations](https://github.com/luksgrin/rust-MAFFT#b-intentional-residuals-tied-trace)
are DP tie-breaks that move a residue to the neighbouring column. For
`--exp 4.30` (FFT-NS-2) on the sample, C and `mafft-rs` give the same
width (517) and 12 `diff` lines (six sequence lines on each side):
513 of 517 columns are identical, 99.9 % of aligned residue pairs are
shared, and MAFFT's `--scoreout` objective is 55 049.75 (C) vs 55 081.77
(`mafft-rs`), a 0.06 % difference. Such a residual is invisible to SP/TC
against a reference unless the moved residue happens to sit in a core
block; byte comparison is what finds it, and what keeps it listed.

### SP / TC of C and Rust on BAliBASE 3

`benchmarks/balibase_scores.py` aligns every BAliBASE 3 set with the
in-tree C MAFFT and with `mafft-rs`, scores both outputs against the
set's XML reference with BAliBASE's own `bali_score` 3.01 (core blocks,
compiled from the archive's `bali_score_src`), and records whether the
two outputs are byte-identical. A Python re-implementation of
`bali_score`'s definitions runs alongside it on every alignment and
agrees with it on all of them. All 386 sets were run in all five modes
(3860 alignments); every C/Rust pair is byte-identical, so every SP and
TC pair is equal. Means over the sets in each subset
(arm64 host, 2026-10-31; raw rows in
`benchmarks/results/balibase_scores_2026-10-31.tsv`):

| Subset | Mode | Sets | Byte-identical | SP C | SP Rust | TC C | TC Rust | SP/TC differ |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| RV11 | FFT-NS-2 | 76 | 76 | 0.577 | 0.577 | 0.308 | 0.308 | 0 |
| RV12 | FFT-NS-2 | 88 | 88 | 0.908 | 0.908 | 0.780 | 0.780 | 0 |
| RV20 | FFT-NS-2 | 82 | 82 | 0.911 | 0.911 | 0.379 | 0.379 | 0 |
| RV30 | FFT-NS-2 | 60 | 60 | 0.852 | 0.852 | 0.500 | 0.500 | 0 |
| RV40 | FFT-NS-2 | 49 | 49 | 0.870 | 0.870 | 0.490 | 0.490 | 0 |
| RV50 | FFT-NS-2 | 31 | 31 | 0.862 | 0.862 | 0.507 | 0.507 | 0 |
| RV11 | FFT-NS-i | 76 | 76 | 0.656 | 0.656 | 0.437 | 0.437 | 0 |
| RV12 | FFT-NS-i | 88 | 88 | 0.926 | 0.926 | 0.828 | 0.828 | 0 |
| RV20 | FFT-NS-i | 82 | 82 | 0.923 | 0.923 | 0.443 | 0.443 | 0 |
| RV30 | FFT-NS-i | 60 | 60 | 0.867 | 0.867 | 0.565 | 0.565 | 0 |
| RV40 | FFT-NS-i | 49 | 49 | 0.890 | 0.890 | 0.542 | 0.542 | 0 |
| RV50 | FFT-NS-i | 31 | 31 | 0.884 | 0.884 | 0.563 | 0.563 | 0 |
| RV11 | L-INS-i | 76 | 76 | 0.679 | 0.679 | 0.463 | 0.463 | 0 |
| RV12 | L-INS-i | 88 | 88 | 0.939 | 0.939 | 0.844 | 0.844 | 0 |
| RV20 | L-INS-i | 82 | 82 | 0.936 | 0.936 | 0.496 | 0.496 | 0 |
| RV30 | L-INS-i | 60 | 60 | 0.871 | 0.871 | 0.607 | 0.607 | 0 |
| RV40 | L-INS-i | 49 | 49 | 0.919 | 0.919 | 0.602 | 0.602 | 0 |
| RV50 | L-INS-i | 31 | 31 | 0.905 | 0.905 | 0.593 | 0.593 | 0 |
| RV11 | G-INS-i | 76 | 76 | 0.706 | 0.706 | 0.482 | 0.482 | 0 |
| RV12 | G-INS-i | 88 | 88 | 0.939 | 0.939 | 0.843 | 0.843 | 0 |
| RV20 | G-INS-i | 82 | 82 | 0.940 | 0.940 | 0.515 | 0.515 | 0 |
| RV30 | G-INS-i | 60 | 60 | 0.879 | 0.879 | 0.619 | 0.619 | 0 |
| RV40 | G-INS-i | 49 | 49 | 0.918 | 0.918 | 0.585 | 0.585 | 0 |
| RV50 | G-INS-i | 31 | 31 | 0.907 | 0.907 | 0.607 | 0.607 | 0 |
| RV11 | E-INS-i | 76 | 76 | 0.667 | 0.667 | 0.435 | 0.435 | 0 |
| RV12 | E-INS-i | 88 | 88 | 0.934 | 0.934 | 0.839 | 0.839 | 0 |
| RV20 | E-INS-i | 82 | 82 | 0.928 | 0.928 | 0.476 | 0.476 | 0 |
| RV30 | E-INS-i | 60 | 60 | 0.866 | 0.866 | 0.594 | 0.594 | 0 |
| RV40 | E-INS-i | 49 | 49 | 0.894 | 0.894 | 0.539 | 0.539 | 0 |
| RV50 | E-INS-i | 31 | 31 | 0.893 | 0.893 | 0.567 | 0.567 | 0 |
| **all** | **FFT-NS-2** | **386** | **386** | **0.826** | **0.826** | **0.500** | **0.500** | **0** |
| **all** | **FFT-NS-i** | **386** | **386** | **0.855** | **0.855** | **0.571** | **0.571** | **0** |
| **all** | **L-INS-i** | **386** | **386** | **0.871** | **0.871** | **0.607** | **0.607** | **0** |
| **all** | **G-INS-i** | **386** | **386** | **0.879** | **0.879** | **0.616** | **0.616** | **0** |
| **all** | **E-INS-i** | **386** | **386** | **0.861** | **0.861** | **0.584** | **0.584** | **0** |

"SP/TC differ" counts sets where C and Rust get a different SP or TC; it
is zero exactly because "Byte-identical" equals "Sets". The absolute
numbers describe MAFFT 7.526, not the port, and are only comparable to
published BAliBASE figures computed the same way (core blocks,
`bali_score` 3.01).

To re-run (alignments are cached under `target/balibase_scores/`, so an
interrupted run resumes; separate invocations merge into the same TSV):

```bash
make -C mafft-upstream/core
cargo build --release -p mafft-rs
export BALIBASE_DIR=/path/to/bb3_release          # RV11/ … RV50/, bali_score_src/
benchmarks/balibase_scores.py run --modes fftns2
benchmarks/balibase_scores.py run --modes fftnsi linsi ginsi einsi --jobs 6
benchmarks/balibase_scores.py summarize benchmarks/results/balibase_scores_<date>.tsv

# the cross-policy comparisons above
CARGO_TARGET_DIR=target/nofma cargo build --release -p mafft-rs --features fp-contract-none
benchmarks/balibase_scores.py crossarch            # BAliBASE, FFT-NS-2
benchmarks/balibase_scores.py counterexamples      # 712/738 and --exp 4.30 on the sample
```

TCS was not computed: T-Coffee is not part of the toolchain, and for
byte-identical outputs it would, like SP and TC, be equal by
construction.

## Byte-identity is per platform build

!!! warning "C MAFFT is not bit-reproducible across compilers"

    This is a property of the original MAFFT, which rust-MAFFT has to
    deal with. Upstream MAFFT 7.526 does not set `-ffp-contract`. The
    same C source therefore produces two alignments:

    - Apple clang fuses `a * b + c` into one rounding. The sentinel
      `--bl 50 --retree 2 --maxiterate 0` on
      `mafft-upstream/test/sample` is 712 columns.
    - gcc, in the ISO C mode the upstream Makefile uses (`-std=c99`,
      so `-ffp-contract=off`), and MSVC, round the multiply and the add
      separately, including on Linux aarch64 where the CPU has FMA
      hardware. The same sentinel is 738 columns.

    A single C binary does not flip between the two. Which one you get
    is fixed by the compiler that built it.

    rust-MAFFT copies that choice when the binary is built, so a default
    install matches the C toolchain used on that platform: fused on
    Apple, not fused on Linux and Windows. The published wheels and
    release binaries are built that way (`*-apple-darwin` fused, Linux
    and Windows wheels not fused). There is no runtime switch, because
    the choice would otherwise be tested on every dynamic-programming
    cell.

    To match the other toolchain, build from source:

    - `cargo build --release -p mafft-rs --features fp-contract-fma`
      matches an Apple-clang C build on any machine
    - `cargo build --release -p mafft-rs --features fp-contract-none`
      matches a gcc or MSVC C build on any machine, including a Mac

    For a from-source Python install, pass the same flag to maturin.
    The steps in
    [Install → From source](../getting-started/install.md#from-source)
    are otherwise unchanged; replace the `maturin develop --release`
    line with one of:

    ```sh
    maturin develop --release --features fp-contract-fma    # Apple clang
    maturin develop --release --features fp-contract-none   # gcc or MSVC
    ```

C MAFFT 7.526 is not bit-identical across CPU architectures. The
upstream Makefile compiles with `-O3` and no `-ffp-contract` setting,
which leaves the floating-point contraction policy to the compiler:

| C reference build | Compiler | Fused multiply-adds per binary |
|--|--|--|
| arm64 (Apple Silicon, `make` in-tree) | clang | on the order of a thousand — 1025 in the census that motivated this section; 882 scalar `fmadd`-family (1139 including vector `fmla`/`fmls`) with Apple clang 21, ~1140 in CI |
| x86-64 (`make` in-tree on Linux) | gcc | 0 |
| arm64 (Apple Silicon, `make CC=gcc-15` in-tree) | Homebrew GCC 15.2 | 0 |

clang fuses `a*b + c` into a single `fmadd` (one rounding); gcc on
x86-64 emits a multiply followed by an add (two roundings). Both are
valid C, and the resulting binaries disagree on some inputs. The policy
follows the compiler, not the architecture: GCC in an ISO mode such as
the Makefile's `-std=c99` does not contract even on arm64, so a GCC build
on Apple Silicon matches the `.nofma` fixtures (checked with
`MAFFT_FP_POLICY=nofma scripts/regen_policy_fixtures.sh … check`) and is
mirrored with `--features fp-contract-none`. The known
sentinel case is `--bl 50 --retree 2 --maxiterate 0` on the 36-sequence
`mafft-upstream/test/sample`:

| Platform | C MAFFT columns | rust-MAFFT columns (default policy) |
|--|--|--|
| arm64 / clang | 712 | 712 |
| x86-64 / gcc | 738 | 738 |

So "byte-identical to C" is always relative to a C build made on the
same platform, and rust-MAFFT tracks that platform rather than picking
one universal answer.

### The reference build

The reference is defined precisely: **the pinned `mafft-upstream` source,
built in-tree with the upstream Makefile's default flags
(`make -C mafft-upstream/core`), on the platform you run on.** No fixture
may be generated from a downloaded binary — a bioconda, Homebrew or
prebuilt-tarball `mafft` was compiled with someone else's toolchain and
flags, and it is not what the tests, the FFI cross-validation or CI compare
against.

### The `fp` policy module

`mafft_types::fp` centralises the decision:

- `mafft_types::fp::CONTRACTS_FMA` — a `const bool`. Defaults to `true`
  on `aarch64` (mirroring clang's behaviour there) and `false`
  elsewhere.
- `mafft_types::fp::fmadd(a, b, c)` — computes `a*b + c` with a single
  rounding when `CONTRACTS_FMA` is set, and as two separate operations
  otherwise. Every hot loop that has to match a contracted C expression
  goes through it rather than calling `f64::mul_add` directly.
- Cargo features `fp-contract-fma` / `fp-contract-none` override the
  default. They are defined on `mafft-types` and forwarded by every
  crate in the workspace, so `cargo test -p mafft-core --features
  fp-contract-none` builds the whole dependency graph under the x86-64
  policy on any host. (`--features` cannot be passed at the
  virtual-workspace root; use `-p <crate>`.)

### Fixture convention

Where a C-canonical fixture differs between the two policies it is
committed twice, as `<name>.fma` and `<name>.nofma`, and the test
helper `fixture_path_fp(name)` picks the one matching `CONTRACTS_FMA`.
Fixtures without a suffix are identical under both policies.

Each pair is described by a line in
`crates/mafft-core/tests/fixtures/policy_fixtures.tsv`
(`fixture_basename<TAB>input_relpath<TAB>mafft_args`) and regenerated by
`scripts/regen_policy_fixtures.sh`, which detects the policy from
`uname -m` (arm64/aarch64 => `.fma`, else `.nofma`; override with
`MAFFT_FP_POLICY=fma|nofma`):

```bash
make -C mafft-upstream/core
scripts/regen_policy_fixtures.sh mafft-upstream/scripts/mafft check   # diff against the checked-in variant
scripts/regen_policy_fixtures.sh mafft-upstream/scripts/mafft write   # overwrite this host's variant
```

A `.fma` fixture therefore comes from the in-tree build on an arm64 host
and a `.nofma` one from the in-tree build on an x86-64 host; never from
the other, and never from a binary you did not build from the pinned
source.

### What CI checks

- The `Build (x86-64)` and `Build (arm64)` jobs each build the C
  submodule in-tree and run the full suite (lib, integration and FFI
  `cross_validate_*` tests) against it, so the claim is checked on both
  reference platforms. Each job prints an FMA census of the C binaries
  it just built and runs the `--bl 50` sentinel above, diffing
  `mafft-rs` against the same-platform C `mafft` wrapper. It then runs
  `scripts/regen_policy_fixtures.sh ... check`, regenerating every
  policy-sensitive fixture from that in-tree C build and failing on any
  byte difference from the checked-in variant for the runner's policy
  (`.nofma` on x86-64, `.fma` on arm64) — so the committed fixtures are
  continuously re-derived from the reference definition rather than
  trusted.
- A `policy-cross-check` job runs the fixture-based (non-FFI) test
  binaries on the arm64 runner with `--features fp-contract-none`, which
  checks the `.nofma` fixtures from an arm64 host. The FFI tests are
  deliberately excluded there: they compare against the C compiled on
  the same runner, which contracts, so under `fp-contract-none` the two
  sides would be running different policies.

## Design decisions that fall out

### No global state

The C code uses ~400 `extern` globals. The Rust modules use owned
`ScoringContext`, `Topology`, `Profile` structs passed explicitly. The
engine is reentrant, embeddable, and trivially thread-safe at the
call-site boundary.

This wasn't done for ergonomic reasons — it was forced by byte-identity.
C MAFFT's globals get assigned in subtle orderings (e.g., `gapfaclocal`
depends on whether `--allowshift` is parsed before or after
`unalignlevel`). Replicating that in Rust required making the
"set-it-now" semantics explicit.

### Hand-ported Cooley-Tukey FFT

`mafft-fft/src/fft_c_compat.rs` is a hand port of MAFFT's `core/fft.c`,
matching the C rounding exactly. Off-the-shelf FFT libraries (rustfft,
realfft, FFTW) vary in butterfly grouping order, producing 1-ULP
correlation differences. On flat-landscape similarity matrices those
ULPs flip FFT-anchor selection — which propagates into different
alignment segmentation and ultimately different output bytes.

The same holds for what is fed to it. For nucleotides C's `Falign`
correlates a single channel (`n20or4or2 = 1`, `Falign.c:923`: only `a`
enters `seq_vec_3`), and the correlation is a count of matching `a`s, so
distinct lags tie exactly and only rounding decides which one `getKouho`
picks first. `mafft_fft::soukan_c_compat` therefore reproduces `Falign`'s
correlation step by step (`calcNaiseki`'s `a · conj(b)` with clang's
contraction, channel sum, inverse FFT, read-back order), and
`cross_validate_falign_dna.rs` checks it bit for bit against C's own `fft`
and `calcNaiseki` on 4–5 kb sequences (issue #43).

### FP order matters

In several hot loops the C compiler's contraction policy decides the
result: with clang on arm64, `a*b + c` becomes a single fused
multiply-add (one rounding); with gcc on x86-64 it stays two operations
(two roundings). rust-MAFFT reproduces whichever the host's C build
does through `mafft_types::fp::fmadd`, which fuses only when
`CONTRACTS_FMA` is set (see [above](#byte-identity-is-per-platform-build)).
One concrete case (`calcW` in `mafft-tree/src/weighting.rs`), matching
clang's `fmadd d2, b, c, d3; fmadd d0, a, b, d2` sequence:

```rust
// a*b + (b*c + a*c), fused exactly where clang fuses it; under the
// x86-64 policy `fmadd` degrades to the unfused mul + add gcc emits.
let s = fmadd(a, b, fmadd(b, c, a * c));
```

Doing this as a naive `a*b + b*c + a*c` on arm64 produces visibly
different weights on `BB40043`, `BB30018`, `BB40010`, `BB30010`,
`BB40004` — and doing it fused on x86-64 would drift from the gcc build
the same way.

### Tie-break ordering

Tied scores in Needleman-Wunsch / Smith-Waterman traceback are resolved
in C via `if (new >= best)` (note the `>=` — equality favors the new
candidate). Several constraint-mode scanners use a different inequality
(`if (new > best)`). Mixing them up flips one trace direction on one
column, which propagates downstream. We mirror each C usage site
exactly.

Sorting ties follow the same rule. `splittbfast` orders sequences with
the C library's `qsort`, whose placement of equal elements is
implementation-defined; `mafft-tree/src/bsd_qsort.rs` is an exact port
of the macOS libc `qsort` (Apple Libc `stdlib/FreeBSD/qsort.c`, with its
`med3`, bounded insertion-sort fallback and `heapsort` depth cap), so the
port matches the macOS C build on every platform. The same applies to the
`rand()` draws in PartTree's pivot selection (macOS `rand()`, unseeded).

### Per-step weight normalization

C MAFFT's `Falign` (the FFT-anchored DP entry) expects the per-group
weight vectors to sum to 1.0. The Rust progressive-alignment caller
normalizes before invocation. Skipping this on a single merge step
produces an off-by-one column on `BB30013`.

## Cross-validation harness

The byte-identity bar is held in place by 82 FFI test functions across
15 files (`crates/*/tests/cross_validate*.rs`; 99 including the three
forensic binaries) that compile MAFFT's C source in-tree via
`mafft-c-bindings` and call both implementations on the same inputs,
comparing outputs byte-for-byte. Because the C is
compiled on the same machine, these tests check Rust's default policy
against the matching C build on whichever platform they run.

This is documented in detail on the
[Cross-validation page](cross-validation.md).

## What this means for changes

Any change to `mafft-core` / `mafft-align` / `mafft-fft` / `mafft-tree`
/ `mafft-scoring` MUST keep the test suite green. New behaviour goes
behind a feature flag or a new function — not as a tweak to an existing
hot path. If you're tempted to "clean up" a slightly weird-looking
arithmetic in the engine, check the cross-validate test for that
function first; it's probably weird because the C output requires it.

## What this DOESN'T mean

- Algorithmic improvements are still possible — they just need their
  own modes / flags, and the byte-identical existing modes must keep
  working.
- Performance optimisation is fine wherever it doesn't change the
  observable output. We've shaved cycles via pre-computed boundary
  tables and LTO+fat codegen — none of which alter the alignment.
  Collapsing arithmetic into `f64::mul_add` is *not* in this category:
  it changes rounding, so it must go through `fp::fmadd` and be
  justified by the C assembly on the matching platform.
- The Python and CLI surfaces can evolve independently of the engine.
  Pretty-printing, progress callbacks, alternate output formats — all
  on the table.
