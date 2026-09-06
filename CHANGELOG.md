# Changelog

All notable changes to rust-MAFFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
from 1.0.0 onwards. Pre-1.0 releases may contain breaking changes between
minor versions; the `mafft` and `mafft-rs` crates aim to keep their
public surfaces stable from 0.1.0 anyway.

## [Unreleased]

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

### Changed

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

[Unreleased]: https://github.com/luksgrin/rust-MAFFT/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.2
[0.1.1]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.1
[0.1.0]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.0
