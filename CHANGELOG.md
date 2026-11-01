# Changelog

All notable changes to rust-MAFFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
from 1.0.0 onwards. Pre-1.0 releases may contain breaking changes between
minor versions; the `mafft` and `mafft-rs` crates aim to keep their
public surfaces stable from 0.1.0 anyway.

## [Unreleased]

### Added

- **Peak-memory and `--thread` benchmarks** (`benchmarks/memory.py`,
  `benchmarks/threads.py`, `benchmarks/pymafft_threads.py`,
  `benchmarks/README.md`; issues #24, #26). `memory.py` records the peak
  RSS of the in-tree C build (the whole `mafft` wrapper pipeline; checked
  against per-engine shims with `--verify-wrapper`) and `mafft-rs` under
  FFT-NS-2, FFT-NS-i, L-INS-i and `--parttree` at `--thread 1/4/8`, plus
  `--memsave`. `threads.py` checks that `mafft-rs` output is byte-identical
  at every thread count and times both sides at `--thread 1/2/4/8/10`;
  `pymafft_threads.py` does the same for `pymafft.align(threads=N)` and
  shows two concurrent calls overlapping (GIL released). Results for an
  Apple M4 on 2026-10-31 / 2026-11-01 in `benchmarks/results/` and in
  `PROFILING.md` ("Memory footprint", "Multi-threading"), the README
  "Performance" section and the Python API page ("Multi-threading from
  Python"). In short: at `--thread 8` `mafft-rs` peaks 2.3–4.8× lower
  than C on FFT-NS-2, FFT-NS-i and L-INS-i (median over 7 inputs per
  mode), but single-threaded it uses about twice C's memory on L-INS-i
  and more on `--parttree`; its output was identical across thread
  counts in every run (incl. all 76 BAliBASE RV11 sets), while its
  `--thread 8` speedup is at most 1.5× (C: 1.6–3.8×) because the
  iterative refinement runs on one thread (#48). `pymafft.align(threads=N)`
  behaves the same (1.05–1.53× at 8 threads, all 90 results equal to the
  CLI at `--thread 1`), and two concurrent calls in two Python threads
  take as long as one (7.53 vs 7.37 s), confirming the GIL is released.
- **High-throughput workload benchmark** (`benchmarks/htworkload_corpus.py`,
  `benchmarks/htworkload.py`, `benchmarks/htworkload/Snakefile`, issue
  #23). A corpus of 2695 *S. pneumoniae* gene families (30 RefSeq
  genomes, MMseqs2 at 70 % identity, 2–100 sequences each) and 1500
  OrthoDB Saccharomycetaceae orthogroups (5–60 sequences), with member
  lists and checksums in `benchmarks/htworkload/`. The driver aligns
  every cluster with `--auto --adjustdirection --thread 1 --nuc` /
  `--auto --thread 1` through Snakemake (C wrapper or `mafft-rs`), a
  plain subprocess loop, and `pymafft.align()` serially and in process /
  thread pools, at 1 / 4 / 8 workers, recording wall-clock, CPU, peak
  RSS and byte-identity with C for the whole run. Pipeline check in
  `benchmarks/results/htworkload_smoke_2026-10-31.tsv`; how to run the
  full workload in `PROFILING.md`.
- **Scaling benchmark corpus and driver** (`benchmarks/scaling_corpus.py`,
  `benchmarks/scaling.py`, issue #22). Ten Pfam families (full-length
  UniProtKB sequences) and RefSeq 16S rRNA, deterministically subsampled
  to n = 100, 250, 500, 1000 (5000 for two families, 2000 for 16S), with
  the manifest, accession lists and checksums in `benchmarks/corpus/`.
  The driver times the in-tree C MAFFT, the `mafft-rs` CLI and native
  `pymafft.align()` interleaved across FFT-NS-2, FFT-NS-i, L-INS-i,
  `--parttree` and G-INS-i at `--thread 1` and 8, recording wall-clock,
  peak RSS, width and byte-identity; per-run caching makes it resumable
  and a projection rule skips cells C cannot finish within an hour.
  Pipeline check in `benchmarks/results/scaling_smoke_2026-10-31.tsv`;
  how to run the full sweep in `PROFILING.md`.
- **Memory-safety comparison** (`benchmarks/memsafety.py`,
  `benchmarks/memsafety_gen.py`, `docs/architecture/memory-safety.md`,
  issue #27). Runs the test corpora, generated edge-case FASTA (long
  headers, long and many sequences, non-ASCII and NUL bytes, empty
  records, line endings, ragged `--add` / `--seed`, malformed `--treein`)
  and 600 mutated copies of the sample through the C reference built with
  `-fsanitize=address,undefined` and through `mafft-rs`; 908 runs in
  `benchmarks/results/memsafety_2026-10-31.tsv`. The C build reports a
  heap-buffer-overflow in `loadtree` for a `--treein` index of 0 or below
  (segfault without sanitizers), a heap-buffer-overflow in `putlocalhom2`
  for `--addfragments --localpair` on valid input, and a null-pointer
  offset in `splittbfast` under `--parttree`; nothing on the other
  inputs. The BAliBASE-wide sweep is deferred; the command is in the doc.

- **BAliBASE SP / TC scores alongside byte-identity**
  (`benchmarks/balibase_scores.py`, issue #28). Aligns every BAliBASE 3
  set with the in-tree C MAFFT and with `mafft-rs` in FFT-NS-2, FFT-NS-i,
  L-INS-i, G-INS-i and E-INS-i, scores both against the core blocks of
  the XML reference with BAliBASE's `bali_score` 3.01 (cross-checked on
  every alignment by a Python re-implementation), and records whether the
  outputs are byte-identical; alignments are cached so runs resume.
  Results in `benchmarks/results/balibase_scores_2026-10-31.tsv`. New
  section "Byte-identity versus alignment scores" in
  `docs/architecture/byte-identity.md`: identical bytes imply identical
  scores but not the converse, shown with the arm64 vs x86-64 C builds
  (712 vs 738 columns on the sample; 8 BAliBASE sets with different bytes
  and the same SP) and the `--exp 4.30` tied-trace residual. The
  `crossarch` and `counterexamples` subcommands reproduce both.
- **Upstream policy** (`docs/architecture/upstream-policy.md`, issue #21).
  Writes down how rust-MAFFT follows upstream MAFFT: the `mafft-upstream`
  submodule gitlink is the reference (currently `0a2319b`, `v7.526` + 2
  build-only commits) and every release names it; the upgrade procedure
  when upstream moves (bump, rebuild the C in-tree, FFI cross-validation
  and full suite, re-derive the C-canonical fixtures with
  `scripts/regen_policy_fixtures.sh` and the fixture READMEs, BAliBASE
  sweep with `scripts/balibase_parity.py`, triage); that fixes to defects
  shared with the C code are reported upstream first; the known residuals
  and why they are carried; and what the port will not do (new default
  behaviour, changed flag semantics, scoring or algorithm changes outside
  upstream). Linked from the architecture index, the home page and the
  contributing guide.

### Changed

- `check-mafft-upstream.yml` now acts on drift instead of only failing the
  run: when the submodule pin is behind upstream HEAD, or upstream has a
  `v*` tag the pin does not contain, it opens an issue titled
  `upstream: mafft-upstream pin is behind <tag-or-describe>` naming the
  pinned and upstream commits/tags, the missing commits and a link to the
  upgrade procedure. An open issue with that prefix is updated (and
  commented on) only when upstream has moved again, so there is at most
  one. The run succeeds once the issue is current. Still weekly plus
  manual dispatch, and it still never moves the submodule; the job now
  needs `issues: write`.
- CI and Python Package workflows now trigger on a positive `paths` list
  (workspace crates, Cargo manifests, the submodule pointer, the scripts
  they call, their own workflow files) instead of a `paths-ignore` list, so
  editing an unrelated workflow or any Markdown no longer starts the
  two-platform build and the wheel matrix. A newer push to the same PR or
  to `main` cancels the run in flight (`concurrency`); release and dispatch
  runs are never cancelled.

### Fixed

- **`--mapout` / `--compactmapout` help text.** Both imply
  `--keeplength` (as in C, `scripts/mafft:699-710`); the help and
  `docs/cli-reference.md` still said they required it. The
  `--addfragments` code comment now describes the `distance=multi` path
  it implements.
- **`--addfragments --reorder` and the linkage options.** `addsingle`
  writes an `order` file (the existing tree's root members, each
  followed by the fragments placed next to it, `addsingle.c:3601-3615`)
  that `f2cl -r` applies; `mafft-rs` ignored `--reorder` on this path.
  `addsingle` also gets the wrapper's `$treealg`, so `--averagelinkage`,
  `--minimumlinkage` and `--mixedlinkage` change both the existing-rows
  tree and the `addonetip` placement (`sueff_global`,
  `mltaln9.c:14464`); `mafft-rs` always used `-X 0.1`.
- **`--anysymbol` restore with `--keeplength`, `--mapout` and
  `--compactmapout`.** C's `restoreu` first blanks the residues
  `--keeplength` dropped from each added sequence (`_deletelist`,
  `restoreu.c:183-206`) and only then copies the original characters
  back; `mafft-rs` copied them in order, so after a dropped insertion
  every residue of that row was shifted. Under `--compactmapout`,
  `restoreu` strips the `Nins:…|` prefix only from the rows that lost
  residues (`insname[i]`, `restoreu.c:270`); `mafft-rs` cut every title
  at its first `|`, so a UniProt title such as `sp|P3|NAME_3` was not
  found and its non-standard residues came out as `X`. Covers `--add`
  and `--addfragments`.
- **`--addfragments --compactmapout` map coordinates.** `addsingle`
  passes the aligned rows to `reconstructdeletemap_compact`
  (`addsingle.c:3819`), so each `vN` is an alignment column, while
  `disttbfast` / `tbfast` (`--add`) pass the raw sequences and `vN`
  counts kept residues. `mafft-rs` used the `--add` convention for both
  (`18v19` where C writes `30v31` on `sample.first36` split 33 + 3). The
  test now compares the `.map` file as well as stdout.
- **`--add --keeplength` keeps the existing alignment's all-gap columns.**
  C records them before the merge and puts them back after the
  insertions are deleted (`recordoriginalgaps` /
  `restoreoriginalgaps`, `disttbfast.c:3909`, `:4862`); `mafft-rs`
  dropped them, so `sample.first36` split 33 + 3 came out 664 columns
  wide instead of 717. This also moves the `--mapout` columns, which are
  counted on the restored rows. `--addfragments -Y` already did this.
- **`--dpparttree` matches C after a profile passes 10,000 columns**
  (issue #45). `splittbfast` `pairalign` switches to `MSalignmm` when
  either profile is longer than 10000 columns, and stays there for the
  rest of that process (`splittbfast.c:1020`). `--parttree` and
  `--dpparttree` share that call (the script does not pass `-F`/`-G`,
  so the post-switch aligner is `MSalignmm`, not
  `Falign_udpari_long`). The Hirschberg port dropped columns when an
  m-split read the midpoint row instead of continuing to `jumpi`, and
  when a split with `jmid` at the end of the second profile skipped the
  half that emits the remaining columns as terminal gaps. A merge is
  now at least as wide as its inputs. PF00069 n = 250 `--dpparttree`
  is 11016 columns and byte-identical to the in-tree C build; PF00400
  n = 250 is 10562 columns and byte-identical as well. `--nomemsave`
  keeps `A__align` past the threshold.
- **The default floating-point contraction followed the ISA, not the C compiler**
  (issue #51). C MAFFT 7.526 leaves `-ffp-contract` unset, so the same
  source produces both alignments: Apple clang fuses `a * b + c` into one
  rounding, and gcc in the Makefile's `-std=c99` mode does not, including
  on Linux aarch64, where the CPU has FMA hardware. MSVC does not contract
  by default either. The 36-sequence `--bl 50 --retree 2 --maxiterate 0`
  sentinel is 712 columns when fused and 738 when not.
  `mafft_types::fp::CONTRACTS_FMA` now defaults to fused only on Apple
  targets, so a Linux aarch64 build (and the published `*-unknown-linux-*`
  wheel) matches gcc instead of Apple clang. `fp-contract-fma` and
  `fp-contract-none` still force either policy at build time; there is no
  runtime switch, and a wheel that should follow the other C build is
  rebuilt from source with the feature. CI on `ubuntu-24.04-arm` builds
  the in-tree gcc reference, checks that it has 0 fused ops, and requires
  the sentinel to be 738 columns and equal to C.
- **Non-UTF-8 FASTA and a CR-only file match C** (part of #41). The reader
  takes the file as bytes. A byte that is not a letter, including one that
  is not valid UTF-8, is dropped and the rest is aligned, as C's residue
  filter does; a sequence name keeps its original bytes (`SeqName`). Every
  `\r` becomes `\n` before records are split, matching the wrapper's
  `tr "\r" "\n"`, so a CR-only file of N sequences yields N sequences. A
  NUL still does not end the line: C keeps only the bytes before it, and
  `mafft-rs` drops the NUL and keeps the rest. That difference is
  intentional.
- The FFI cross-validation tests did not compile on Linux aarch64: 56
  sites in 15 test files spelled the C `char**` argument (and the `alg`
  global) as `i8`, which is `c_char` on x86-64 Linux and Apple arm64 but
  not on Linux aarch64, where `char` is unsigned. They now use `c_char`. CI
  builds and tests the workspace on `ubuntu-24.04-arm` on each push.
- **FFT modes diverged from C on long nucleotide sequences** (issue #43).
  C's `Falign` correlates nucleotides on a single channel — `n20or4or2 =
  1` when `scoremtx == -1` (`Falign.c:923`), so only `a` enters the FFT —
  while the port correlated all ten scored letters. That ranks the 20
  candidate lags differently; on similar sequences the best lag carries
  the anchors either way, but on three 3.6–5.0 kb *S. pneumoniae* CDS
  (`dna_01715`) the second progressive merge picked different anchors
  under FFT-NS-1/2/i and `--auto`. The FFT anchor search now also follows
  C's arithmetic exactly: the correlation is C's `calcNaiseki` (`a ·
  conj(b)`) read back in `Falign`'s order (`mafft_fft::soukan_c_compat`,
  bit-identical to C on this input), the per-site anchor scores use
  `alignableReagion`'s summation order and its division by `totaleff`,
  the FFT length is the power of two strictly above the longer profile
  (as `Falign` rounds), and the twiddle table recurrence is fused like
  arm64 clang's (`make_sintbl`; identical below 2^19 points). Not
  sensitive to the FP contraction policy (C built with
  `-ffp-contract=off` gives the same bytes). New FFI tests
  `cross_validate_falign_dna.rs` (correlation bit for bit, candidate
  lags, and every merge of the cluster replayed through C's `Falign`)
  and the `dna_01715` end-to-end fixture.
- **`--add --reorder`, `--add --maxiterate` and `--compactmapout` names match C** (issue #50). `--reorder` follows `topolorderz` of disttbfast's second guide tree (`-E 2`), not the k-tuple tree that drove the merge. `--maxiterate N` runs `dndpre -M 2` and `dvtditr -K nadd -B` after the progressive add, with the single-threaded rule when `--thread` is absent. `--compactmapout` implies `--keeplength` and prefixes added names that lost residues with `Nins:<first-run>|`.
- **`--addfragments` and `--addfragments --keeplength` match C** (issue #50). The default path is `addsingle` on `distance=multi` (`scripts/mafft:2628`): `dndpre -M 2` distances among the existing rows, `pairlocalalign -Y` for the fragment distances and local homology, then one constrained profile merge (`-l 2.7`, `importance = opt * 0.5`, new gaps `=`). `--keeplength` is `addsingle -Y`. On the 30 + 6 sample the widths are 871 and 595. `--addfragments --localpair` stays unpinned: C reads off the end of a buffer on some of those inputs.
- `mafft-rs` lost residues without an error when a profile alignment had
  one empty side: an input record with no residues printed every
  sequence empty, and DNA with a long homopolymer (e.g. 300 × `A`) lost
  residues under FFT-NS-1/2/i and `--memsave`. The empty side is now
  placed against gaps, as C `A__align` does (issue #27).
- `--add` / `--addfragments` onto rows of different lengths panicked
  (index out of bounds); it now exits 1 with C's "The original N
  sequences must be aligned" message.
- A malformed `--treein` (index 0, negative, out of range, too few lines)
  ended the process from inside the engine with `process::exit(1)`, also
  for in-process callers; it is now a normal `MafftError`.
- **`--add` printed `=` in the alignment** (issue #40). When a new-gap
  region borders restored common gaps, C's `insertnewgaps` realigns it
  with `profilealignment`, which marks the existing side with `=`; C
  clears those marks after every `--add` merge (`eq2dashmatometehayaku`,
  `disttbfast.c:2988`, `tbfast.c:1705`). The port skipped that step, so
  later merges aligned `=` as a residue and it reached the output. It is
  now cleared as in C; the 30 + 6 reproducer and its 2 + 1 minimisation
  are byte-identical to C, as are `--keeplength` and `--anysymbol` on
  the same input. `mafft-rs` also checks the final alignment and exits 1
  with an internal-error message, instead of printing it, if a row holds
  `=`, `<` or `>` (bytes no input residue can be).
- **`--parttree` panic and divergence on 50+ sequences** (issue #38).
  The PartTree port covered only the first level of
  `splittbfast.c::splitseq_mq`. From 50 sequences on (`picksize` no
  longer exceeds `njob`, so `tokyoripara = 0.7`) redundant pivots are
  dropped and each yuko group holds distinct sequences that C aligns in a
  recursive `splitseq_mq` call before merging; the Rust side merged them
  unaligned, and `merge_step_cached` indexed past the shorter rows. New
  `mafft-tree::parttree_splitseq` ports the full recursion (pivot draw
  with the macOS libc `rand()` sequence, redundancy filter, duplicate
  handling, `treeorder` recursion, `pairalign` merge order and
  orientation) for the k-tuple pass, the `-Z` (`fromaln`) pass and the
  `-U` DP pass. `--parttree --reorder` / `--treeout` and `--dpparttree`
  use it; `--dpparttree` now runs one pass with DP distances, as C's
  `DPPartTree-1` does. PartTree merges no longer use FFT (`splittbfast`
  is never given `-F`). `bsd_qsort` is now an exact port of the macOS
  libc `qsort` (its `med3` tie rule, bounded insertion-sort fallback and
  `heapsort` depth cap differed for tied elements), and the
  `splittbfast` length factor follows the FP contraction policy.
  Byte-identical to the in-tree C build on a 50-sequence regression
  fixture, 82 of 84 scaling-corpus cells (10 Pfam families × n = 100 and
  250 × `--parttree` / `--reorder` / `--nofft` / `--dpparttree`, 16S
  `--nuc`); the two `--dpparttree` residuals at n = 250 are listed in
  `docs/architecture/byte-identity.md`.

- The in-tree C reference (`mafft-c-bindings/build.rs`) builds with GCC 14
  and 15 (issue #17). GCC 14 made `-Wincompatible-pointer-types` an error
  that `-w` does not silence, and `iteration.c` passes `char (*)[256]`
  to `writePre(char **, ...)` twice, so `cargo test` failed in the build
  script. `build.rs` now passes `-Wno-error=` for that diagnostic and for
  the other GCC 14 / Clang 16 promotions (implicit function declaration,
  implicit int, int conversion, return mismatch, incompatible function
  pointer types), each only if the compiler accepts it. With Homebrew
  GCC on macOS, `build.rs` also links GCC's `libemutls_w.a`, which the
  shim needs for `__thread` and which rustc's clang link step does not
  provide. `make -C mafft-upstream/core` was not affected (it pins
  `-std=c99` and never compiles `iteration.c`). CI gains a
  `C build: gcc 14` / `gcc 15` job that builds the shim and runs the
  upstream Makefile in the official `gcc` images.

### Performance

- **`--thread N` scales the iterative modes** (issue #48). The per-call
  rayon pool is unchanged. On it, independent units write into
  preallocated slots and are consumed in a fixed order, so `mafft-rs` at
  every thread count matches `mafft-rs --thread 1`, which matches C
  `--thread 1` (C's `athread` rules when `nthread > 0`). C's own
  `--thread N` for N ≥ 2 is nondeterministic and is not the reference.
  No DP cell arithmetic and no decision-boundary summation changed
  (`mafft_types::fp::fmadd` operand order is the same). What runs on the
  pool:
  - the all-vs-all pairwise alignments of L-INS-i, G-INS-i and E-INS-i
    (already parallel; unchanged here);
  - progressive merges of disjoint guide-tree subtrees, one level at a
    time (`--c-compat` stays serial, because its cpmx memo carries from
    merge to merge), and the FFT segments of one profile-profile merge;
  - the FFT-NS-i refinement segments, each started from its own slice;
  - refinement branches: up to `--thread N` branches are scored
    speculatively against the same alignment, then applied in the
    `--thread 1` visiting order. An accepted branch discards the later
    speculative results, which are re-scored against the updated
    alignment. `--thread 1` still selects C's `athread` convergence rule
    and does no speculative work.
  Two single-thread changes sit on the same path: `build_imp_matrix`
  maps a raw column to its residue index in constant time, and the
  `G__align11` row loop hoists the substitution-matrix row and skips the
  warp-record copies when warp is off.
  Idle-machine wall-clock medians for `--thread 1` versus `--thread 8`
  are not in this entry yet; they replace this sentence after the formal
  run.

### Python

- `pymafft` builds against PyO3 0.30 (was 0.24), which supports
  Python 3.14. The deprecated `Python::with_gil` / `allow_threads` calls
  are now `Python::attach` / `detach`. The build-system floor is
  `maturin>=1.9.4`, and Python 3.14 is in the classifiers (#18).
- CI: `python.yml` gains a `lint` job that runs
  `RUSTFLAGS="-D warnings" cargo check -p pymafft`, so new PyO3
  deprecations fail CI instead of piling up. Publishing waits for it.
- `pymafft` wheels target the stable ABI (PyO3 `abi3-py39`): one
  `cp39-abi3` wheel per target installs on CPython 3.9 and every later
  version, so a release ships 5 wheels plus the sdist instead of 25, and a
  new Python minor no longer needs a rebuild. `requires-python` stays
  `>=3.9`. `python.yml` builds each target once and installs that same
  wheel on Python 3.12 (push / pull request) or 3.9, 3.12 and 3.14
  (release / dispatch) (#19).

## [0.2.0] - 2026-09-06

**Output changes since 0.1.x.** Every item below moves rust-MAFFT *toward*
the C MAFFT 7.526 reference, but each one changes bytes or exit codes for
anyone comparing against 0.1.x output:

- Nucleotide alignments are written in lowercase (C's read-time case fold);
  protein stays uppercase; `--anysymbol` / `--preservecase` still keep the
  input's own case.
- DNA pairwise gap penalties use C's `3 *` nucleotide scale, so L-INS-i /
  G-INS-i / E-INS-i and `--auto` on nucleotide input can align differently.
- DNA refinement guide trees use the nucleotide `dndpre` offset (matrix
  shift 220, not the protein 73), so FFT-NS-i and every refinement mode
  without a pairwise phase can change on DNA.
- Exactly-two-sequence inputs are now refined once, unweighted, under
  `--maxiterate N > 0`, every `*-INS-i` mode and `--auto` (they were
  returned unrefined before).
- `--thread N` with N >= 1 follows C's `athread` refinement rules
  (ascending branch walk, per-cycle convergence, `Converged2.` /
  `Oscillating?`), so refined output under `--thread N >= 1` changes; no
  `--thread` and `--thread 0` are untouched.
- Illegal residues are rejected with `Illegal character c` and exit 1
  (`pymafft.MafftError`, code 1) as C's `seqcheck` does: protein `U` / `O`,
  nucleotide `.`, anything outside the alphabet, in the main input and the
  `--add` file. A blank before `>` is rejected with C's format message and
  exit 1. `.` is kept as a legal protein residue instead of being stripped.
- Type detection counts `N` as nucleotide, so N-rich DNA is no longer
  aligned as protein.
- `searchAnchors` flushes the trailing anchor at the loop's exit index, so
  FFT-NS-i's last refinement segment no longer starts one column early when
  the final conserved region runs to the end of the alignment.
- Floating-point contraction now follows the reference C build of the
  target: fused on aarch64, not fused elsewhere. x86-64 builds therefore
  match x86-64 C MAFFT for the first time (`--bl 50 --retree 2 --maxiterate 0`
  on the 36-sequence sample gives width 738 on x86-64 and 712 on arm64; 0.1.x
  gave 712 everywhere). `--features fp-contract-fma` / `fp-contract-none`
  override the default.

The `--nuc` / `--amino` flags, the `run_from` / `MafftError` / `Mafft` /
`Progress` library API, the nucleotide case fold, the DNA pair-phase gap
scale, two-sequence refinement, the nucleotide `dndpre` offset, the `athread`
convergence rule, `MAFFT_RS_REFINE_STATS` and the per-alphabet constant audit
below were contributed by Johan Henriksson (@mahogny) and integrated from
mahogny/rust-MAFFT for issue #1, with the original authorship preserved.

### Added

- **Policy-sensitive fixture regeneration** (`scripts/regen_policy_fixtures.sh`,
  `crates/mafft-core/tests/fixtures/policy_fixtures.tsv`). The reference for
  every `<name>.fma` / `<name>.nofma` fixture is now defined as the pinned
  `mafft-upstream` source built in-tree with the upstream Makefile's default
  flags on the platform you run on — never a downloaded binary. The script
  regenerates each manifest entry from that build and either diffs (`check`)
  or overwrites (`write`) the variant for the host's policy (`uname -m`,
  overridable with `MAFFT_FP_POLICY=fma|nofma`); CI runs `check` on both the
  x86-64 and the arm64 runner after the `--bl 50` sentinel. The upstream
  drift check (`check-mafft-upstream.yml`) now also runs weekly and still
  only reports.
- **Floating-point contraction policy** (`mafft_types::fp`). C MAFFT 7.526
  is not bit-reproducible across CPU architectures: the arm64 macOS binaries
  contain ~1030 `fmadd`/`fmsub` instructions each (clang contracts `a*b+c`
  into a single-rounding FMA, `scripts/fma_census.sh` counts them), while
  baseline x86-64 builds (gcc `-O3`, same in-tree source) contain none because baseline
  x86-64 has no FMA unit. On the 36-seq protein sample with `--bl 50 --retree
  2 --maxiterate 0` the arm64 C binary gives width 712 and the x86-64 C
  binary gives 738. rust-MAFFT previously hard-coded `f64::mul_add` (always
  fused, and emulated in slow software on x86-64), so it matched only the
  arm64 build; @mahogny's fork removed every `mul_add` and matched only the
  x86-64 build. All 74 library call sites now go through
  `mafft_types::fp::fmadd(a, b, c)`, which is `a.mul_add(b, c)` when
  `mafft_types::fp::CONTRACTS_FMA` is `true` and `a * b + c` otherwise, with
  operand order and nesting preserved from the C source. The policy is
  chosen by cargo features, forwarded by every workspace crate:
  `fp-contract-fma` (always fuse), `fp-contract-none` (never fuse), both is
  a compile error, neither means "mirror the reference C build of the
  target you run on": fused on `aarch64`, not fused elsewhere. Sites that
  arm64 clang deliberately does *not* contract (e.g. `sequence_weights`'s
  `rootnode[s] += len * eff[s]`) stay plain on every target. Fixtures whose
  bytes depend on the policy are committed as `<name>.fma` (arm64 C) and
  `<name>.nofma` (x86-64 C) with the plain name removed; `fixture_path_fp`
  in the mafft-core tests picks the variant matching `CONTRACTS_FMA`. Only
  `sample.bl50.fftns2` needed the split (the `.nofma` copy is the x86-64
  fixture from mahogny/rust-MAFFT). The FFI test that compares the step-24
  BL50 profile DP against the in-tree C build skips under
  `fp-contract-none`, since the C it compiles on an arm64 host is fused.
- `--nuc` / `--amino`: force the input sequence type, overriding the
  ATGC-frequency auto-detection (matches C MAFFT `scripts/mafft:547-550`,
  `seqtype="-D"` / `seqtype="-P"`). Mutually exclusive; inert when absent,
  so auto-detection is unchanged for every existing command line.
- `mafft_rs::run_from(argv, out)`: the argv-driven, non-exiting form of
  `run()`. Same clap definition and therefore exactly the same flag
  semantics (`--auto`'s size heuristic, `--adjustdirection`'s strand
  detection, …), but every path where `run()` calls `std::process::exit(N)`
  returns a `MafftError` carrying the same message text and exit code, and
  the alignment is written to a caller-supplied `io::Write` instead of
  stdout. `--output FILE` still writes to that file.
- `mafft_rs::MafftError`, with `code()` and `message()`.
- `mafft_rs::Mafft`: a typed builder that constructs an argv and hands it
  to `run_from`, so it cannot drift from the command line.
- `mafft_io::apply_case_convention`: applies C MAFFT's residue-case fold
  (lowercase nucleotide, uppercase protein) to a parsed `SequenceSet`.
- `MAFFT_RS_REFINE_STATS=1` prints one line per iterative-refinement call to
  stderr: `refine: nseq=.. len=.. cycles=n/max visited=.. changed=..
  accepted=.. exit=maxiter|converged|oscillation`, plus a
  `refine-segments: anchors=.. segments=..` line for the segmented
  (FFT-NS-i) path. C's `dvtditr` reports its refinement work directly
  (`Segment n/N`, then one line per branch), so this makes "did both sides
  run the same cycles?" answerable without a debugger — a speed comparison
  is meaningless otherwise. Off by default; CLI output is unchanged.
- Progress sink: `mafft_rs::Progress` (one method, `message(&self, &str)`)
  with `StderrProgress` (current behaviour) and `SilentProgress`, plus a
  blanket impl so any `Fn(&str)` is a sink.
  `mafft_rs::run_from_with_progress(argv, out, &(dyn Progress + Sync))`
  routes the run's 12 progress messages there instead of stderr, and
  `Mafft::progress(sink)` does the same for the builder. Both default to
  stderr, so `run_from`, `run` and the CLI are unchanged. The sink is taken
  by shared reference and called through `&self`, so one sink can serve a
  whole worker pool. Only progress is routed — failures come back as
  `MafftError`, and non-fatal `Warning:` / `Could not …` diagnostics plus
  `--scoreout`'s score line stay on stderr, so a silent sink cannot hide
  either a problem or requested output. The sink covers every progress
  message `mafft-rs` emits; two further lines in `mafft-core` (the
  Q-INS-i/X-INS-i BPP line and the `--skipiterate` banner) are not routed —
  see the `progress` module docs for why.
- `mafft_rs::run_from_seqs(argv, &SequenceSet, &(dyn Progress + Sync)) ->
  Result<MultipleAlignment, MafftError>` and `Mafft::run_seqs(&SequenceSet)`:
  the in-memory form of `run_from`. Same clap definition, same flag layer
  (`--auto`'s size heuristic, `--adjustdirection`, `--nuc` / `--amino`,
  `--reorder`, `--thread`, `--add FILE`, …), but the primary input is a
  borrowed `SequenceSet` instead of the positional INPUT file and the
  alignment comes back as a `MultipleAlignment` (names + aligned rows by
  value, plus guide tree / distance matrix / step trace) instead of FASTA
  text. The rows are byte-for-byte what `run_from` prints for the same
  sequences and flags — same order, `_R_` names, case and gap characters —
  which `tests/inmemory_parity.rs` asserts on the 36-sequence protein sample,
  the 120 x 1400 bp CDS fixture and synthetic reverse-complement /
  mixed-case inputs. Input is normalised the way the FASTA reader would
  (residue filter, case convention, type detection when `seq_type` is
  `Unknown`) and copied only when that changes something; a canonical set
  reaches the engine borrowed. `--output FILE` is still honoured; the
  `--treeout` / `--distout` side files need an INPUT path and warn as they do
  for stdin. Errors are the same `MafftError`s as `run_from`. Internally
  `run_from_with_progress` is now `parse_argv` → `preflight` → `read_input`
  → `align_prepared` → `write_alignment`, with `run_from_seqs` swapping
  `read_input` for `prepare_in_memory`; the CLI's bytes are unchanged
  (verified with `cmp` against the pre-change binary and arm64 C MAFFT on
  eight command lines). `mafft-rs` re-exports `MultipleAlignment`,
  `Sequence`, `SequenceSet` and `SeqType` so a caller needs no other crate.
- `mafft` gains an opt-in `cli` feature that re-exports the `mafft-rs` crate
  as `mafft::cli`, so a library user can reach `run_from` / `run_from_seqs`
  / `Mafft` from the umbrella crate. Off by default because it brings clap.
- `mafft_io::normalize_residues`, `residues_are_normalized` and
  `residues_follow_case_convention` expose the FASTA reader's own residue
  filter and case rule for residues that never went through a file, and say
  whether applying them would change anything; `normalize_residues` returns
  the case-preserving reader's `= < >` error so an in-memory caller fails
  where a file would. `mafft_io::find_illegal_residue` is C `seqcheck`'s
  alphabet test (`locaminod` / `locaminon`, taken from `mafft-scoring`, on
  which `mafft-io` now depends); `mafft-rs` runs it in `align_prepared`, so
  a file and an in-memory `SequenceSet` are refused with the same
  `Illegal character c` (exit 1). `detect_seq_type` (and
  `detect_seq_type_with_limit`, now exported) accept any iterator of byte
  slices, so the readers no longer copy every row into a `Vec<Vec<u8>>`
  just to detect the type. Existing `&Vec<Vec<u8>>` callers are unchanged.

- **DNA parity fixtures** with expected outputs captured from the in-tree
  reference build: the five panaroo-rs gene clusters
  (`crates/mafft-bin/tests/fixtures/panaroo_tiny_dna_clusters/`, run as
  `--auto --adjustdirection --thread 1 --nuc` and without `--thread`), the
  120 × 1.4 kb CDS reproducer (`mtb_cds_120x1400.fa`, `.thread1` and
  `.fftnsi` variants), its first 8 sequences (`mtb_cds_first8.fa`), and the
  minimal reproducers `dna_pair_gapscale_min.fa` / `refine_njob2_min.fa`.
  Consumed by `crates/mafft-bin/tests/c_parity_dna.rs`.
- **Input-handling corpus**: 22 inputs × {default, `--anysymbol`,
  `--localpair --maxiterate 0`, `--nuc`, `--add` modes} = 71 cells in
  `crates/mafft-bin/tests/fixtures/input_handling/`, with C's stdout as
  `.expected` and C's exit-1 cases pinned by status + message
  (`tests/input_handling_vs_c.rs`; the evaluation notes are in that
  directory's `README.md`).

### Changed

- CI: the build/test job is the reusable `build-test.yml`, run once per
  reference platform (ubuntu x86-64, macOS arm64). `python.yml` builds one
  wheel per target for Python 3.12 on push / pull request and the full
  5 targets × Python 3.9–3.13 matrix on release / dispatch. `ci.yml`,
  `python.yml` and `docs.yml` run on pull requests against *any* base
  branch; `ci.yml` and `python.yml` skip documentation-only changes
  (`docs/**`, `**/*.md`, `mkdocs.yml`, `CITATION.cff`, `LICENSE*`).
- `--thread N` now builds a *local* rayon pool and `install()`s the
  alignment into it instead of calling `build_global()`. A process-global
  pool can only be initialised once, so an in-process caller running many
  alignments was previously stuck with the first call's thread count.
  No change for the CLI.
- `run()` is now a thin wrapper around `run_from` (unchanged signature and
  observable behaviour: same stdout, stderr, exit codes and messages).
- The FASTA path no longer clones every aligned row into the output
  `SequenceSet`; the rows are moved out of the `MultipleAlignment`. Output
  bytes are unchanged.

### Performance

- **Flattened the profile-DP and `L__align11` scratch buffers.** In
  `profile_align_imp_multimtx` (mafft-align `profile.rs`) the sparse
  `cpmxpd`/`cpmxpdn` views of both profiles are now one flat entry list
  plus a `len + 1` offsets table each instead of `Vec<Vec<(usize, f64)>>`,
  the `h`/`ijp` DP matrices are flat row-major `(n + 1) * (m + 1)` buffers
  with stride `m + 1` (thread-local pools of `Vec<f64>`/`Vec<i32>`, grow-only
  as before), `match_calc_row_into` walks the scoring matrix j-outer /
  l-inner over contiguous rows and consumes the cpmx2 entries sequentially
  through exact-length slices and iterators, and the boundary-init `bcarr`
  buffer is pooled instead of allocated twice per call. In `local_align`
  (`local.rs`) a thread-local `LocalScratch` pool replaces seven allocations
  per call, `seq2` is mapped to alphabet indices once, the match score is
  added inside the DP cell (same operands, same order) instead of in a
  separate row-fill pass, and the j-loop body is the const-generic
  `local_dp_row::<HAS_ROW>` kernel over `m + 1`-cell row slices. No
  `unsafe` was added and every multiply-add still goes through
  `mafft_types::fp::fmadd` with unchanged operand order, so output is
  byte-identical to C MAFFT 7.526 (arm64) on the 36-sequence protein sample
  (default, `--maxiterate 1000`, `--localpair`/`--globalpair`/`--genafpair
  --maxiterate 1000`, `--bl 50 --retree 2 --maxiterate 0`, `--parttree`) and
  on the DNA reproducers (`mtb_cds_120x1400.fa --retree 2 --maxiterate 2`,
  1742 columns), and the `fp-contract-none` build still reproduces the
  x86-64 width 738 on `--bl 50`. Wall-clock medians of 5 interleaved runs on
  Apple Silicon (M4): `mtb_cds_120x1400 --retree 2 --maxiterate 2` 4.99 s ->
  4.62 s (1.08x), sample `--localpair --maxiterate 1000` 0.85 s -> 0.77 s
  (1.11x), sample `--maxiterate 1000` 0.54 s -> 0.48 s (1.14x), sample
  default 0.066 s -> 0.055 s (1.21x). Extracted from the data-layout /
  allocation hunks of mahogny/rust-MAFFT 0b3955d ("fixes and optimization")
  by Johan Henriksson (@mahogny) and re-applied on top of the `fmadd`
  migration.

### Fixed

- **The input read path now matches C MAFFT's.** Four rules were wrong or
  missing, found by running the 22-file corpus above against C (23 of 71
  cells differed; 1 remains — see Notes): (1) type detection now counts `N`
  as nucleotide (`io.c` `countATGC`), so N-rich DNA is no longer aligned as
  protein; (2) C's `seqcheck` (`mltaln9.c:60-85`) is applied after reading
  and to the `--add` file, so protein `U` / `O` and any letter outside the
  alphabet exit 1 with `Illegal character c` (banner unless `--quiet`)
  instead of being aligned silently, and nucleotide `.` is fatal as in C;
  (3) `.` is never stripped — it is a legal protein residue (index 23) and
  C's `gappick0` / `commongappick` remove `-` only, which also fixes the
  L-INS/G-INS/E-INS pair phase seeing gapped input; (4) `--anysymbol`
  follows `charfilter` (drop only `\n`, space, `\r`; `= < >` inside a
  sequence exit 1), `replaceu` (digits / tabs kept as `X` / `n` residues)
  and `restoreu` (originals copied onto every non-`-` position), including
  for the `--add` file. A blank before `>` is rejected with C's format
  message and exit 1 instead of gluing the header into the previous
  sequence. `mafft_io::IoError` gained two variants carrying C's messages.
- **FFT-NS-i's last refinement segment could start one column early.**
  `searchAnchors` (C `mltaln9.c`) slides a 20-column window over the
  alignment with `for( i=1; i<len-divWinSize; i++ )` and, when a
  high-conservation region is still open at the end, closes it with
  `seg->end = i` — the loop's *exit* value `len - divWinSize`. The port
  closed it with the last *iterated* column (`len - divWinSize - 1`), which
  puts the trailing anchor one column to the left whenever `start + len` is
  even. Segments 1..n-1 were unaffected, so this only showed on inputs whose
  final anchor region runs to the end of the alignment: the first 8
  sequences of `mtb_cds_120x1400.fa` under `--retree 2 --maxiterate 1000`
  (C anchor 1412, ours 1411; one gap column moved in `s5`, with and without
  `--nofft` / `--thread 1`), while 5-7, 9-12, 20, 30 and all 120 sequences
  were already identical. Pinned by
  `c_parity_dna::fftnsi_first8_dna_trailing_anchor_matches_c`
  (`crates/mafft-bin/tests/fixtures/mtb_cds_first8.fa`, expected bytes from
  the in-tree arm64 build) and a unit test on the flush index.
- **DNA pairwise gap penalties were a third of C MAFFT's, so L-INS-i /
  G-INS-i / E-INS-i diverged on nucleotide input.** C scales the pair-phase
  gap penalties by `3 * 600/1000` for nucleotide and `600/1000` for protein
  (`constants.c:316-322` vs `:672-677`), keeping the offset at `1 * 600/1000`;
  the pair phase here applied the protein factor unconditionally (a variable
  named `scale_protein`). Gaps were too cheap, so the local/global pairwise
  step bought extra matches with gaps C refuses, and the `hat3` constraints
  and final alignment followed. Now `pair_penalty_scales(is_nucleotide)`,
  with a unit test pinning the `3 *`. Protein and FFT-NS-2 were never
  affected (protein has no `3 *` in C either; FFT-NS-2 does not use these
  penalties).

  **DNA pairwise alignments will change.** Any `--localpair`, `--globalpair`,
  `--genafpair` or `--auto` run on nucleotide input can now produce a
  different alignment than before. As with the case fold, this is a
  correction *toward* the C MAFFT 7.526 reference the crate claims
  byte-identity with, not a behaviour of our own: on 60 synthetic clusters
  L-INS-1 went from 24/60 to 60/60 byte-identical with C, and on 30
  clusters evolved from real biological ancestors `--auto` went from 14/30
  to 30/30. Minimal reproducer: two 15 bp sequences under
  `--localpair --maxiterate 0` (`crates/mafft-bin/tests/fixtures/dna_pair_gapscale_min.fa`).
- Two latent instances of the per-alphabet-constant class, found by the
  audit recorded under Notes rather than by a parity failure. Neither was
  reachable from the engine, so no output changes: `GapModel::default()` was
  protein-shaped *and* arithmetically wrong (`-918`; C truncates toward zero,
  giving `-917`) and is now correct and documented as protein-only; and the
  **public** `FftAlignParams::dna()` inherited that protein gap default
  instead of DNA's `-2753`, which would have mis-scaled gaps for any
  downstream caller using it. `GapPenalties::default()` holds unscaled
  `ppenalty`-style units unlike every other `GapPenalties` in the tree; it is
  unused, and now says so.
- **`--thread N` (N >= 1) now follows C's `athread` refinement rules.** C
  picks its refinement implementation on `nthread > 0` (`tditeration.c:1433`),
  and with one worker `athread` is deterministic but differs from the
  single-threaded `TreeDependentIteration` in four ways, all now modelled
  behind `RefinementParams::per_cycle_convergence` (selected by
  `MafftEngine::nthread`, matching C's `-C` mapping: no `--thread` and
  `--thread 0` give `-C 0`, the single-threaded rules):
  * *Branch order.* The single-threaded loop reverses the walk on odd cycles
    (`:1641-1648`); `athread` consumes `branchtable[0..nbranch]`, the identity
    permutation under the default `randomseed = 0` (`:522`, `:728-730`), so
    every cycle walks the tree in ascending order. This was the last 2-line
    `--thread 1` residual on the 120-sequence reproducer: one gap column in
    `s97` moved because cycle 1 was walked backwards.
  * *Convergence.* Single-threaded tests `converged >= locnjob * 2` after
    **every branch** and stops mid-cycle (`:2328-2342`); `athread`'s collector
    tests once per **cycle** whether any branch gained (`maxgain > 0.0`,
    `:590`) and only stops at the top of the next cycle (`:527-551`). C's own
    output shows it: at `--maxiterate 2`, 22 of 85 segments print
    `Converged.` alone, 56 print `Converged.` *and* `Reached 2`, 7 print
    `Reached 2` alone.
  * *Oscillation.* Single-threaded compares a branch's score with the same
    branch 2, 4, 6 … cycles earlier and exits at once (`:2345-2372`).
    `athread` instead raises `Converged2.` when a branch's score equals its
    score in **any** cycle `<= iterate-2` (`:1217-1230`, `:636-637`), and
    `Oscillating?` when the cycle's last accepted score equals an earlier
    cycle's (`:609-619`); both stop at the end of the cycle. On the 36-seq
    protein sample `--thread 1 --localpair --maxiterate 1000` C stops via
    `Converged2.` in cycle 4 and rust-MAFFT now does the same.
  * *Skipped branches* (`--skipiterate`) still record `tscore = mscore`
    (`:1064-1067`, `:1234`) so they take part in the `Converged2.` check.

  **DNA output changes for `--thread N >= 1` with refinement.** The
  120-sequence FFT-NS-i reproducer is byte-identical to C `--thread 1` at
  `--maxiterate 1/2/3/1000`, with and without `--retree 2`, and to C's
  `--thread 2` output on this input; the 36-seq protein sample is
  byte-identical under `--thread 1` in FFT-NS-i, L-INS-i, G-INS-i and
  E-INS-i at `--maxiterate 1000`; the synthetic cluster corpus under
  `--auto --adjustdirection --thread 1 --nuc` is 60/60. The no-`--thread`
  path is untouched and remains byte-identical to C. Pinned by
  `fftnsi_120seq_dna_matches_c_under_thread_1`
  (`mtb_cds_120x1400.thread1.expected`, captured from the in-tree arm64
  build).
- **DNA refinement guide trees used the protein `dndpre` offset, reordering
  UPGMA merges.** For modes with no `pairlocalalign` step (FFT-NS-i and
  friends) the refinement tree is rebuilt the way C's `dndpre` does. C's
  script does not pass `-h` to that `dndpre` call, so `constants()` uses its
  DEFAULT `poffset` — and the alphabets do not share one: `DEFAULTOFS_N =
  -369` (`DNA.h:3`) gives a matrix shift of 220, `DEFAULTOFS_B = -123`
  (`blosum.c:3`) gives 73. The shift was hardcoded to the protein 73.

  On DNA that produced refinement distances differing from C's `hat2`
  outright (on a 120-sequence 1.4 kb input, leaf pair (0,58): C 0.253, ours
  0.307), which swapped two UPGMA merges, which changed the group on 131 of
  237 refinement branches and so the final alignment. Now
  `dndpre_offset_shift(is_nucleotide)`, with a unit test pinning both values
  against the C constants.

  **DNA output changes for FFT-NS-i and any refinement mode without a
  pairwise phase.** Byte-parity with C MAFFT 7.526 over BAliBASE `bali2dna`
  (141 real DNA benchmark sets) under `--maxiterate 2` goes **67/141 → 132/141**,
  and the 120-sequence reproducer becomes byte-identical. Protein is
  unaffected (it already used 73), as are FFT-NS-2 and `--localpair`, which
  never reach this path.
- **Two-sequence inputs were never refined.** Every refinement entry point
  returned early at `nseq <= 2`. C does not skip a pair: `dvtditr.c:704-708`
  sets `weight = 0; niter = 1` for `njob == 2`, `tditeration.c:772` then
  uses uniform weights, and `:1425` gates branch-weight computation on
  `locnjob > 2`. So a pair is refined exactly once, unweighted — which can
  change it (e.g. `AA`/`CC` under `--maxiterate 1000`: C gives `aa`/`cc`,
  we gave `aa-`/`-cc`). Guards relaxed to `nseq < 2` and the iteration cap
  mirrors `niter = 1`; `BranchWeights` already yielded uniform weights at 2.
  Affects `--maxiterate N > 0` and every `*-INS-i` mode, including `--auto`,
  on exactly-two-sequence input, DNA and protein alike.
- **Nucleotide output case now matches C MAFFT.** DNA/RNA alignments were
  emitted in uppercase; C MAFFT emits them in lowercase. C folds residue
  case as it reads — `io.c:1462-1467` (`load1SeqWithoutName_realloc`) calls
  `onlyAlpha_lower` when `dorp == 'd'` and `onlyAlpha_upper` otherwise, and
  `readData_pointer` repeats the nucleotide pass with `seqLower`
  (`io.c:1755`; its `upperCase != -1` guard is only reachable from the
  legacy non-FASTA `FRead` header parser, so it is always true for FASTA
  input). rust-MAFFT's reader applied `onlyAlpha_upper` unconditionally.
  It now folds per the sequence type, and `--nuc` / `--amino` re-apply the
  fold for the type they force, because C's `$seqtype` fixes `dorp` before
  any sequence is read.

  **This changes existing output**: a DNA/RNA alignment that previously
  came back uppercase now comes back lowercase. That is the point — it is
  what C MAFFT 7.526 produces, and whole-file case flips were breaking
  byte-for-byte comparisons against a C-MAFFT reference. Protein output is
  unchanged (uppercase, as before and as in C). `--anysymbol` /
  `--preservecase` are unchanged too: they keep the input's own case, in
  both C and here.

  This closes the parity caveat the README recorded as
  "byte-exact (case-insensitive)"; `--nofft samplerna` is now byte-exact
  with `cmp`, not just with `diff -i`.
- `--adjustdirection` / `--adjustdirectionaccurately` help text claimed the
  k-mer strand detection was "not yet implemented"; it has been implemented
  since TODO R-5 (2026-06-03).

### Python

- `pymafft` now runs every alignment through the CLI's flag layer
  (`mafft_rs::run_from_seqs`) instead of driving `MafftEngine` directly. Each
  keyword is one `mafft-rs` flag, so nothing about a flag's meaning is
  re-derived in Python and the output is byte-identical to the command
  line's for the same flags (82 option-parity tests compare `to_fasta()`
  against the binary's stdout on the same fixtures). `align`, `align_file`
  and `align_fasta_string` keep `strategy` / `maxiterate` with their
  previous defaults and gain keyword-only options: `strategy="auto"`
  (`--auto`), `seq_type="nuc"|"amino"` (`--nuc` / `--amino`; `None` detects
  as the CLI does), `adjust_direction=True|"accurately"`
  (`--adjustdirection` / `--adjustdirectionaccurately`), `threads`
  (`--thread`), `reorder`, `retree`, `scoring="bl62"|"jtt200"|"tm100"`
  (`--bl` / `--jtt` / `--tm`), `gap_open` / `gap_extend` (`--op` / `--ep`),
  `quiet` (default `True`) and `progress` (a callable receiving the lines
  the CLI prints on stderr; with `quiet=False` and no callable they go to
  stderr as before). The GIL is released while the alignment runs.
- `pymafft.run(args, sequences, *, progress=None)`: escape hatch that
  passes an explicit flag list straight to `run_from_seqs`, for flags
  without a keyword.
- `pymafft.MafftError`, a `ValueError` subclass carrying the CLI's exit code
  (`.code`) and stderr text (`.message`); `except ValueError` keeps working.
- `AlignmentResult.to_biopython()` returns a `Bio.Align.MultipleSeqAlignment`;
  Biopython is imported lazily, only inside that call.
- **Two behaviour changes inherited from the reader fixes above**, because
  Python input now takes the same path as a FASTA file on the command line:
  (1) nucleotide alignments come back lowercase (they were uppercase) and
  protein input is uppercased on read, so `align(["acdefghik", …])` returns
  `ACDEFGHIK`; (2) residues outside the alphabet (protein `U` / `O`,
  nucleotide `.`, …) raise `MafftError` (`Illegal character U`, code 1)
  instead of being aligned as unknowns — as C MAFFT does; use
  `run(["--anysymbol", …], seqs)` to keep them. Two tests were updated for
  this (`test_mixed_lengths`, and `test_case_passes_through_verbatim`
  became `test_case_follows_c_mafft`).
- `pymafft` depends on the `mafft-rs` crate (path dependency); its
  `fp-contract-*` features are forwarded there too.

### Notes

- **Per-alphabet constant audit.** Three bugs were found reactively where a
  constant correct for one alphabet was applied on a path serving both
  (`scale_protein` in the pair phase, the `dndpre` offset, and the two latent
  ones above). Every value in the translation deriving from C's
  `constants()`, `DNA.h`, `blosum.c` or `JTT.c` has now been enumerated and
  checked against *both* C branches (nucleotide `constants.c:296-326`,
  protein `:664-682` / `:895-910`) — 17 live sites, all correct.
  `mafft_scoring::penalties` carries a test asserting every one of them for
  both alphabets at once, including that the `offset*` values take C's `1 *`
  factor on the nucleotide branch while the gap penalties take `3 *`, so a
  future edit cannot give one alphabet the other's constant unnoticed.

- C MAFFT 7.526 genuinely produces different output for *no* `--thread` than
  for `--thread 1` — it selects a different refinement implementation on
  `nthread > 0` (`tditeration.c:1433`), and the two converge by different
  rules and walks the tree in a different order. This was previously recorded
  here as "C-side sensitivity" that rust-MAFFT could not match; that was
  wrong, and rust-MAFFT now reproduces both paths byte-for-byte (see the
  `--thread` entry under Fixed). Both are deterministic: 5/5 identical over
  repeat runs.
- Still genuinely unmatchable: `--thread N` for **N >= 2**. C is
  nondeterministic there — the same binary on the same input produced 2
  distinct outputs over 3 runs at both `--thread 2` and `--thread 4` — so
  byte-identity with C is impossible in principle at those thread counts.
  rust-MAFFT remains deterministic across all thread counts.
- One input-handling cell still differs from C: with `--add --anysymbol`,
  an added sequence ending in `*` (→ `X`) is placed `…KWRR-----X` by C and
  `…KWRRX-----` here. It reproduces with a literal trailing `X` and no
  `--anysymbol`, so it is a tie-break in the `--add` profile DP on a
  zero-scoring terminal residue, not an input-handling bug; the test
  `add_anysymbol_gapped_existing_with_unusual_new` is `#[ignore]`d with
  that reason (residual R-A in `TODO.md`).
- **FFI cross-validation tests serialised; poison-tolerant guards.** The
  test binaries that link the C reference (`mafft_sys`) call functions that
  read/write ~150 process-wide C globals (and `treeCnv` keeps a static scratch
  buffer), while cargo runs `#[test]` functions on multiple threads. Six FFI
  test files (`cross_validate_weights`, `cross_validate_counteff`,
  `cross_validate_cpmx`, `cross_validate_bb20027_dp`,
  `exp_residual_dump_compare`, `cross_validate_insertnewgaps`) had no
  serialisation and could crash intermittently (`cross_validate_weights`
  SIGABRT/SIGSEGV on macOS arm64 CI); they now take the same `C_MUTEX` as the
  other twelve. Every guard is now `unwrap_or_else(PoisonError::into_inner)`
  so one genuine failure no longer cascades into misleading `PoisonError`s.
  The scoring binary's `init_c_globals()` also resets `TMorJTT = JTT`:
  `initglobalvariables()` does not reset it, so a `*_tm_*` test running first
  leaked the TM table into the C-default JTT probe (the intermittent "JTT
  mismatch"). Test-only changes; no crate `src/` touched.

## [0.1.2] - 2026-06-10

Metadata fixes. No engine, library, or CLI behaviour changes from
0.1.1 / 0.1.0; bumped so that all distribution channels can publish
together cleanly.

### Fixed

- `CITATION.cff`: dropped the SPDX expression `MIT AND BSD-3-Clause`
  in favour of the single SPDX identifier `MIT` (Zenodo's `cffconvert`
  pipeline rejected the expression form with "Citation metadata load
  failed"). The BSD-3-Clause attribution for the algorithmic constructs
  ported from upstream MAFFT remains in `LICENSE-BSD` and the workspace
  `Cargo.toml` `license = "MIT AND BSD-3-Clause"` (where cargo accepts
  expressions fine).
- Maintainer email updated to `lucas.goiriz@csic.es` in workspace
  authors and `pymafft` pyproject metadata.

### Notes

- v0.1.1 published partially: docs + GH-release binaries succeeded;
  crates.io got 5 of 10 crates before hitting the new-crate rate
  limit; PyPI publish was gated by a transient quay.io docker-pull
  flake on one Linux wheel. v0.1.2 retries all channels with a valid
  CITATION.cff.

## [0.1.1] - 2026-06-10

Release-pipeline fixes. No engine, library, or CLI behaviour changes
from `0.1.0`; bumped only because `pymafft 0.1.0` was published to
PyPI before the rest of the release-day workflows could be fixed.

### Fixed

- `release.yml`: added `permissions: contents: write` so
  `softprops/action-gh-release@v2` can attach the cross-compiled
  binaries to the GitHub release (previously failed with
  "Resource not accessible by integration").
- `docs.yml`: the `github-pages` deployment environment is now scoped
  to allow tag-triggered deploys (`v*`) in addition to `main` pushes.
- `python.yml`: x86_64-apple-darwin wheel is now cross-compiled from
  the Apple Silicon `macos-latest` runner (the previously-targeted
  `macos-13` Intel runner pool was retired by GitHub in 2026).
- `python.yml`: Linux wheel build no longer changes the cwd inside the
  manylinux container (`cd ../..` had broken maturin's manifest
  resolution; replaced with `--manifest-path ../../Cargo.toml`).

### Notes

- `pymafft 0.1.0` was yanked on PyPI after this release. Use 0.1.1+
  for any new installs. Programmatic and CLI behaviour is unchanged.

## [0.1.0] - 2026-06-10 [YANKED]

Initial release attempt. `pymafft` published to PyPI successfully but
the crates.io publish (email verification not yet completed), GitHub
release binary upload (missing `contents: write` permission), and
docs deploy (environment-rule rejected the tag) all failed. Superseded
by 0.1.1.

### Engine

- **Byte-identical to C MAFFT 7.526** across the full BAliBASE 3
  fixture set (1930/1930) for every supported mode: FFT-NS-2, FFT-NS-i,
  G-INS-i, L-INS-i, E-INS-i, PartTree, DPPartTree, plus `--add`,
  `--adjustdirection` (k-mer mode), `--oneiteration`, `--bestfirst`,
  `--allowshift`, `--unalignlevel`, `--skipiterate`, `--pileup`,
  `--youngestlinkage` / `--averagelinkage` / `--minimumlinkage` /
  `--mixedlinkage`.
- **Stub-parity** for `--pdbidlist` and `--pdbfilelist`: print the
  upstream "temporarily unavailable, 2018/Dec." message and exit 0 to
  match C MAFFT 7.526 behaviour verbatim (these were disabled
  upstream).
- **Pure Rust runtime**: the release binary compiles zero C code.
  `mafft-c-bindings` (the FFI cross-validation shim) is a dev-only
  internal crate.

### Distribution

Four parallel release channels off a single GitHub release tag:

- **Rust library**: `cargo add mafft` (top-level re-export of
  `mafft-core`, `mafft-types`, `mafft-io`). Sub-crates also published
  individually for finer-grained dependence.
- **Rust CLI**: `cargo install mafft-rs`.
- **Python package**: `pip install pymafft` ships 25 wheels (Linux
  x86_64 + aarch64, macOS Intel + Apple Silicon, Windows x86_64 ×
  Python 3.9–3.13) plus an sdist.
- **Python executable**: the `pymafft` wheel bundles the `mafft-rs`
  binary inside it; `pip install pymafft` puts `mafft-rs` on `$PATH`
  via a console-script entry. No Rust toolchain required.
- **Pre-built CLI binary**: attached to every GitHub release for 5
  targets (Linux x86_64 + aarch64, macOS Intel + Apple Silicon,
  Windows x86_64).

### Testing

- ~430 Rust tests including ~140 FFI cross-validation tests that
  compile MAFFT C source in-tree and compare per-function outputs
  byte-for-byte.
- 76 Python tests across binding parity, biopython interop, and the
  bundled-CLI smoke path.
- Full CI matrix on every push / PR: Linux byte-identity gate, plus
  macOS + Windows build + lib-test cross-platform sanity.

### Documentation

- Site at https://luksgrin.github.io/rust-MAFFT (Material for MkDocs)
- Per-crate `cargo doc` published to https://docs.rs

[Unreleased]: https://github.com/luksgrin/rust-MAFFT/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.2.0
[0.1.2]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.2
[0.1.1]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.1
[0.1.0]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.0
