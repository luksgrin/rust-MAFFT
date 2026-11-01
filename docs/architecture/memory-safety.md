# Memory safety and robustness

This page records what happened when the project's test corpora, and a set
of malformed or unusual FASTA inputs of the kind users pass by accident,
were run through two builds: the C reference compiled with
AddressSanitizer and UndefinedBehaviorSanitizer, and `mafft-rs` (issue #27).
It lists each sanitizer report, what `mafft-rs` does on the same input, and
the inputs on which nothing was reported.

The subject is **undefined behaviour and robustness**: out-of-bounds
access, crashes, panics, silent data loss. It is not memory *footprint*.
Peak RSS and allocation volume are issue #26.

## Method

**C side.** The pinned `mafft-upstream` source, built in-tree with clang
and the standard sanitizer flags, then installed into `target/asan/`:

```bash
nice -n 19 make -C mafft-upstream/core -j6 CC=clang STRIP=true \
    CFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' \
    LDFLAGS='-fsanitize=address,undefined' \
    PREFIX=$PWD/target/asan/mafft-bin install
make -C mafft-upstream/core clean
```

Each run goes through upstream's own `mafft` wrapper script with
`MAFFT_BINARIES=target/asan/mafft-bin/libexec/mafft`,
`ASAN_OPTIONS=detect_leaks=0:abort_on_error=0:log_path=…` and
`UBSAN_OPTIONS=print_stacktrace=1:log_path=…`. The `log_path` setting
matters: under `--quiet` the wrapper sends engine stderr to `/dev/null`,
so a report printed to stderr would be lost. Leak detection is off because
a leak at exit is not undefined behaviour. ASan stops a process at its
first error. UBSan reports and continues.

**Rust side.** `cargo build --release -p mafft-rs`, same arguments, same
120 s timeout.

**Harness.** `benchmarks/memsafety.py` renices itself to 19. For each run
it records the C exit status, stdout size (the wrapper can exit 0 with
empty stdout after an engine error), sanitizer report counts, the first
report line with its first `mafft-upstream` frame, shell-reported signals
and a stderr excerpt. On the Rust side it records exit status, outcome
(`ok`, `error`, `panic`, `signal`, `timeout`) and a stderr excerpt. When
both succeed it also records byte-identity and whether the two outputs
hold the same residues per record (`residues`). That last check separates
gap-placement differences from lost or altered data. Results are in
[`benchmarks/results/memsafety_2026-10-31.tsv`](https://github.com/luksgrin/rust-MAFFT/blob/main/benchmarks/results/memsafety_2026-10-31.tsv)
(908 runs, arm64 macOS, Apple clang 21). Re-running resumes.

```bash
benchmarks/memsafety_gen.py          # generated inputs → target/memsafety/
benchmarks/memsafety.py              # all small corpora
benchmarks/memsafety.py --summary
```

## Corpora

| corpus | runs | content |
|---|---:|---|
| `input_handling` | 110 | the 22 files of `crates/mafft-bin/tests/fixtures/input_handling/` × default / `--anysymbol` / `--localpair --maxiterate 0` / `--nuc` / `--add add_new_unusual.fa`, including the inputs where C exits 1 by design |
| `sample` | 22 | `mafft-upstream/test/sample` × 15 modes (FFT-NS-1/2/i, G/L/E-INS-i, `--parttree`, `--dpparttree`, `--memsave`, `--adjustdirection[accurately]`, `--auto`, `--thread 2`, `--treeout`, `--anysymbol`); `samplerna` × 7 |
| `clusters` | 72 | the 5 panaroo DNA clusters, the 30 `r2_dna_clusters` and the 120 × 1.4 kb CDS fixture × `--auto --adjustdirection --thread 1 --nuc` and FFT-NS-i |
| `sample_add` | 7 | 30 rows of the C sample alignment + the last 6 sequences: `--add`, `--addfragments`, `--keeplength`, `--reorder`, `--maxiterate 2`, `--localpair` |
| `adversarial` | 94 | `memsafety_gen.py`, seed 27: headers of 255–1,000,000 bytes, a 120 kb sequence (wrapped and on one line), a 30 kaa protein, 3000 short and 2000 tiny sequences, UTF-8 / Latin-1 / UTF-16 / random bytes, NUL bytes, empty records (first, middle, last, all), empty file, no header, CR-only / CRLF / mixed line endings, missing final newline, duplicate names, gap-only / `*`-only / digit-only sequences, `--add` onto ragged rows (± `--keeplength`, `--addfragments`), empty / gap-only / over-long / non-ASCII added sets, `--seed` ragged / single / empty / gap-only, 14 malformed `--treein` files, and a long low-complexity DNA set for `--adjustdirection` |
| `mutants` | 600 | `mafft-upstream/test/sample` with byte substitutions, truncations, inserted `>` / `-` / NUL, deleted newlines, duplicated headers, `>` at line start (75 per kind, seed 29), FFT-NS-2 only |
| `cases` | 3 | `benchmarks/memsafety_cases/`: minimised inputs from the sweep |

An earlier pass over 2000 mutants of the same generator (before the
`mafft-align` fix below) also produced no sanitizer report. The final table
keeps 600 so the whole sweep stays within about 15 minutes of niced CPU.

## Findings: sanitizer reports from the C build

| # | input | C (ASan / UBSan) | C without sanitizers | `mafft-rs` |
|---|---|---|---|---|
| C1 | `--treein` file with an index of `0` or a negative index | heap-buffer-overflow, READ of size 4 in `loadtree` `mltaln9.c:2636` | `Segmentation fault: 11` in `disttbfast`, exit 1 | exit 1: `--treein: line 1: im=0 or jm=2 out of range [1, 6]` / `invalid im '-5'` |
| C2 | `--addfragments … --localpair` on a valid input (1 existing sequence + 2 fragments, minimised from the sample) | heap-buffer-overflow, READ of size 4 in `putlocalhom2` `io.c:803` (called from `pairlocalalign`) | exit 0, deterministic output (3/3 runs) | exit 0, every residue kept, output differs from C |
| C3 | `--parttree` / `--dpparttree` on `sample` and `samplerna` | `splittbfast.c:2377:101: runtime error: applying non-zero offset 8 to null pointer` | exit 0 | exit 0, byte-identical to C |

**C1.** `loadtreeoneline` (`mltaln9.c:1482`) decrements both indices and
rejects only `ar[0] >= ar[1]`. `loadtree` then rejects `im > nseq-1` but
not `im < 0`, so a `0` (or negative) index reaches `hist[im]`:

```
==98970==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x603000001e6c
READ of size 4 at 0x603000001e6c thread T0
    #0 in loadtree mltaln9.c:2636
    #1 in disttbfast disttbfast.c:3953
0x603000001e6c is located 4 bytes before 24-byte region
allocated by thread T0 here:
    #1 in AllocateIntVec mtxutl.c:239
    #2 in loadtree mltaln9.c:2525
```

The same read happens under FFT-NS-i (`--maxiterate 2`). The parser in
`mafft-tree/src/treein.rs` checks both bounds. Until this change it was
called from inside the engine, which reported the error with
`std::process::exit(1)`. That also ended in-process callers such as
pymafft. `mafft-rs` now validates the tree before aligning and returns a
normal error. Regression tests: `treein_zero_index_is_rejected` and
`treein_negative_index_is_rejected` in
`crates/mafft-bin/tests/memsafety_regressions.rs`.

**C2.** A legitimate command line: `--localpair` is one of the pair modes
the wrapper accepts with `--addfragments`:

```
==16948==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x60b0000005bc
READ of size 4 at 0x60b0000005bc thread T0
    #0 in putlocalhom2 io.c:803
    #1 in pairalign pairlocalalign.c
    #2 in pairlocalalign pairlocalalign.c:3258
0x60b0000005bc is located 4 bytes before 104-byte region
allocated by thread T0 here:
    #1 in AllocateIntMtx mtxutl.c:300
    #2 in constants constants.c:889
```

The 104-byte region is one row of `n_dis` (26 ints). The read is
`n_dis[amino_n[*pt1]][amino_n[*pt2]]` with the second index `-1`. That fits
a byte outside the alphabet, most likely the terminating NUL of the shorter
of the two strings walked in step. The fault needs at least two fragments.
It was reduced from the 30 + 6 sample set to one existing sequence and two
fragments (`benchmarks/memsafety_cases/addfragments_localpair_*.fa`). It
does not occur with `--add` instead of `--addfragments`. Without
sanitizers, C reads whatever sits before that row and finishes. Its output
is stable across runs but depends on that read, so it is not a reference.
`mafft-rs` places two fragment residues differently. Regression test
`addfragments_localpair_completes_and_keeps_residues` pins that `mafft-rs`
completes, keeps every residue and produces equal-width rows. It does not
pin bytes.

**C3.** `splitseq_mq( …, children+ii, … )` with `children == NULL`.
Pointer arithmetic on a null pointer is undefined. In practice the pointer
is not dereferenced on this path and the output matches `mafft-rs` byte for
byte.

Per the [bug-fix direction](upstream-policy.md#bug-fix-direction), C1–C3
were reported upstream with these reproducers on 2026-10-31:
C1 → [sysimm/mafft#8](https://gitlab.com/sysimm/mafft/-/issues/8),
C2 → [sysimm/mafft#6](https://gitlab.com/sysimm/mafft/-/issues/6),
C3 → [sysimm/mafft#7](https://gitlab.com/sysimm/mafft/-/issues/7).

## Findings: `mafft-rs`

The sweep found more defects on the Rust side than on the C side. None is
undefined behaviour. Three were silent wrong results.

| # | input | before | now |
|---|---|---|---|
| R1 | an input record with no residues (empty record, gap-only, `*`-only, digits-only; a duplicated header line) | **every sequence printed empty**, exit 0 | byte-identical to C (gap-only row for the empty record) |
| R2 | DNA with a long homopolymer, e.g. 300 × `A` + 300 nt, FFT-NS-1/2/i or `--memsave` | **residues silently dropped**, exit 0 (11 of 1195 in the fixture; 1254 of 29,476 in the low-complexity set) | every residue kept; gap placement inside the run still differs from C (608 vs 600 columns) |
| R3 | `--add` / `--addfragments` onto rows of different lengths | panic: index out of bounds, `add.rs` `commongappick` | exit 1 with C's message, `The original N sequences must be aligned` |
| R4 | `--add` of the last 6 sample sequences onto 30 aligned rows | **`=` characters in the output** (3 rows, 5 each), exit 0 | byte-identical to C (fixed in #49) |

**R1 and R2** share one cause. `profile_align`, `fft_profile_align` and
`msalignmm` returned no operations when one profile had zero columns, so
the other profile's columns disappeared from the merge. An empty record
reaches that case directly. A homopolymer reaches it through FFT anchors
inside the run, which leave anchor segments with one empty side. C's
`A__align` fills the empty side with gaps and returns 0
(`Salignmm.c:1228-1255`). `Alignment::against_gaps` now does the same. On
the 600 mutants this turned the 75 duplicated-header inputs from "all
residues lost" into byte-identical output. Tests:
`empty_record_matches_c` and `homopolymer_dna_keeps_every_residue`. A third
test, `homopolymer_dna_matches_c`, is `#[ignore]`d because the remaining
gap-placement difference is open. The full workspace suite passes with
the fix (541 passed, 0 failed).

**R3.** C checks the row lengths (`disttbfast.c:4147-4160`,
`addsingle.c:3217-3228`). `addsingle` exits 1. `disttbfast` jumps to
`chudan`, after which the wrapper exits **0 with empty stdout**. `mafft-rs`
prints C's message and exits 1 on both paths. Stdout is empty, as in C.
Only the exit status differs from the `disttbfast` path. Test:
`add_to_unaligned_existing_is_rejected`.

**R4** is fixed in #49 (issue #40). When a new-gap region borders
restored common gaps, C's `insertnewgaps` realigns it with
`profilealignment`, which marks the existing side with `=`. C clears
the marks after every `--add` merge (`eq2dashmatometehayaku`,
`disttbfast.c:2988`, `tbfast.c:1705`). The port in
`mafft-core/src/progressive.rs` skipped that step, so later merges
aligned `=` as a residue and it was printed. It is now cleared as in C.
`mafft-rs` also refuses to write an alignment holding `=`, `<` or `>`,
none of which the reader admits as a residue: it exits 1 with an
internal-error message. The reproducer is
`benchmarks/memsafety_cases/add_sentinel_*.fa`, and its `--add`,
`--add --keeplength` and `--add --anysymbol` output is now
byte-identical to C. It minimises to 2 existing rows + 1 new sequence.
Tests: `add_sentinel_matches_c` (30 + 6) and
`add_sentinel_minimal_byte_identical_to_c` (2 + 1).

## No diagnostic

ASan and UBSan reported nothing for any other run:

- all 110 `input_handling` runs. C and `mafft-rs` agree on success or
  failure in every cell, and all 66 successful cells are byte-identical;
- all 72 cluster runs, all byte-identical, including the 120 × 1.4 kb set
  under FFT-NS-i;
- 19 of the 22 `sample` / `samplerna` runs (the other 3 are C3), all
  byte-identical;
- the 600 mutants, and the earlier 2000;
- in the adversarial set: headers up to 1 MB, the 120 kb sequence, 3000 and
  2000 sequences, NUL and non-ASCII bytes, binary garbage, UTF-16, empty
  and header-only files, every line-ending variant, duplicate names,
  gap-only rows, `--add` with ragged / empty / gap-only / over-long /
  non-ASCII sets, all `--seed` variants, and the 11 malformed `--treein`
  files other than the 0 / negative index.

No sanitized C run was killed by a signal, and no run on either side timed
out. After the R3 fix no `mafft-rs` run panicked.

## Other differences seen (not memory safety)

These are documented here so they are not lost. None is undefined
behaviour. The NUL case is an improvement over C. The rest are
robustness or parity gaps.

!!! note "A zero byte in the file is kept"

    This is an improvement, and it stays. C reads the record into a C
    string and filters it with `onlyAlpha_*` and `strcpy`
    (`io.c:1355-1377`, `io.c:1751`), which stop at the first zero byte,
    so every residue after it in that record is silently discarded,
    including those on later lines. `AC\0DE` is aligned as `AC`. `mafft-rs`
    drops only the zero byte, because it is not a letter, and aligns
    the residues after it, so the same input is `ACDE`. A zero byte in
    a name is kept, together with the bytes that follow it; C cuts the
    name at that byte (`myfgets`, `io.c:1082`). A normal FASTA file contains no zero byte, so the two
    programs still match on real input. Of the mutants whose residues
    differed in the 2026-10-31 sweep, 74 of 76 were this case.

- **Non-UTF-8 input.** The FASTA reader takes the file as bytes. A
  non-alphabetic byte, including one that is not valid UTF-8, is dropped
  and the remaining residues are aligned, which is what C does. A name
  keeps its original bytes. The 2026-10-31 run refused these files
  (`stream did not contain valid UTF-8`): 57 of 75 substitution mutants
  and every non-ASCII and binary input. That refusal is gone.
  `input_bytes_vs_c` compares a Latin-1 byte inside a sequence, and a
  non-UTF-8 name, with the in-tree C build.
- **Line endings.** Every `\r` becomes `\n` before records are split,
  matching the wrapper's `tr "\r" "\n"`. CRLF and a CR-only file of N
  sequences both yield N sequences. The CR-only case is compared with
  the in-tree C build in `input_bytes_vs_c`. Previously a lone CR was
  not a line break, so that file read as one record.
- **Header length.** C keeps about the first 250 bytes of a name.
  `mafft-rs` keeps all of it, so long names give different bytes but the
  same alignment.
- **Validation C does not do, and the reverse.** C accepts a `--treein`
  that merges the same node twice (exit 0, no diagnostic), which
  `mafft-rs` rejects. `mafft-rs` accepts a `--seed` with ragged rows or a
  single sequence, and a `--treein` line longer than C's 999-byte
  `fgets` buffer, all of which C rejects. C's wrapper also exits 0 with
  empty stdout for an empty file, a file without a header line and NUL- or
  garbage-only input.
- **Byte-identity outside the parity corpora.** Both succeed but the bytes
  differ for: the 3000 × 30 aa set (width 121 vs 129), the 120 kb sequence
  (a few gap moves), the low-complexity `--adjustdirection` set, and `--treein`
  with `nan`/`inf` branch lengths. Of the `sample_add` runs, `--add` and
  `--add --keeplength` have matched since #49. `--add --maxiterate 2`,
  `--add --reorder` and the `--compactmapout` added-sequence names match
  since #50, on the canonical 30 + 6 fixture and on the #40 reproducer.
  `--addfragments` and `--addfragments --keeplength` match as well
  (`addsingle` on the `distance=multi` branch): widths 871 and 595 on
  the 30 + 6 fixture. `--addfragments --localpair` is still not a
  reference (C2 above).

## Published reports

Checked on 2026-10-31:

- [GitLab `sysimm/mafft` issues](https://gitlab.com/sysimm/mafft/-/issues)
  has four issues in total. None is an input-triggered crash. #3, output
  nondeterminism with `--auto --thread 10`, is threading rather than
  memory safety, and its input is behind a login-only upload. #4 is a
  GCC 15 build failure, already handled by this repository.
- [GSLBiotech/mafft#9](https://github.com/GSLBiotech/mafft/issues/9) /
  [steineggerlab/ufcg#42](https://github.com/steineggerlab/ufcg/issues/42)
  report a `tbfast` segfault in 7.526 L-INS-i on 687 sequences ×
  4532 aa. The input is not published, so it cannot be reproduced.
- Biostars ([297087](https://www.biostars.org/p/297087/),
  [344646](https://www.biostars.org/p/344646/)) and the
  [QIIME 2 forum](https://forum.qiime2.org/t/problem-running-mafft-alignment/6193)
  report allocation failures and segfaults on 7.2xx–7.4xx with tens of
  thousands of sequences. Those are memory footprint (#26), not
  reproducible on 7.526 as reported.
- The [MAFFT changelog](https://mafft.cbrc.jp/alignment/software/changelog.html)
  lists for 7.511 "`--adjustdirection` could abort when input sequences
  have long low-complexity region". An input of that kind is in the
  adversarial set (`adjustdirection_lowcomplexity`). C 7.526 does not abort
  and the sanitizers report nothing. It is the input that exposed R2.

## Not covered

- **BAliBASE.** The 386-set sweep under the sanitizers was deferred so it
  would not compete with timing benchmarks. Command, three modes, resumable:

  ```bash
  nice -n 19 benchmarks/memsafety.py --corpus balibase \
      --balibase "$BALIBASE_DIR" --jobs 4 \
      --out benchmarks/results/memsafety_balibase_$(date +%F).tsv
  ```

- **One platform.** arm64 macOS only. Linux / glibc and x86-64 builds may
  lay out the heap differently, which matters for what an out-of-bounds
  read returns.
- **What the sanitizers cannot see.** Uninitialised reads (MemorySanitizer
  is not available on macOS), overflows inside one allocation (C's
  fixed-size `char name[B]` arrays inside larger blocks) and data races
  (ThreadSanitizer was not run; `--thread N ≥ 2` appears only once).
- **Fuzzing depth.** The mutants are one bounded deterministic pass, not
  coverage-guided fuzzing. The Rust side has `cargo fuzz` targets. A
  coverage-guided campaign against the sanitized C engines would go
  further.
