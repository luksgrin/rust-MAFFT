# Input handling vs C MAFFT 7.526 — evaluation of fork hunks B-1 and B-2

Branch `issue1/input-hunks` (base `origin/main` = 3e54a54). Reference: arm64 C
MAFFT `v7.526 (2024/Apr/26)` at `/Users/apophis/.local/bin/mafft`. Source
citations are into the `mafft-upstream` submodule (`core/`), file:line.

The two hunks under study come from fork commit 0b3955d,
`crates/mafft-bin/src/lib.rs`:

* **B-1** — a `thread_local! OnceCell<Option<rayon::ThreadPool>>` cached
  1-thread pool used for `--thread 1` instead of a `ThreadPoolBuilder` per call.
* **B-2** — `replace_unusual` applied to *every* input unconditionally (main
  applies it only under `--anysymbol`/`--preservecase`) plus a new
  `strip_input_gaps` that deletes `-` **and `.`** from input sequences before
  alignment (skipped under `--anysymbol` or `--add`).

## 1. What C does on the read path

### 1a. Default (no `--anysymbol`)

| Step | Where | Behaviour |
|---|---|---|
| Line endings | `mafft.tmpl:1132` `cat "$1" \| tr "\r" "\n" > infile` | CRLF → LF before any binary sees the file. |
| Header format check | `mafft.tmpl:1827-1834` | `grep -c '^[[:blank:]]\+>'` > 0 → prints "The first character of a description line must be / the greater-than (>) symbol, not a blank. / Please check the format around the following line(s):" + `grep -n` and **exit 1**. Printed even under `--quiet`. |
| Type detection | `io.c:2518-2570` `getnumlen` → `io.c:2065-2090` `countATGC` | Among `isalpha` chars, `a t g c u` **and `n`** count as nucleotide; `> 0.75` → `dorp='d'`. Gaps and non-alpha are not counted. |
| Residue filter | `io.c:1435-1470` `load1SeqWithoutName_realloc` → `io.c:1355-1378` `onlyAlpha_lower/upper` | Keep only `isalpha(c) || c=='-' || c=='*' || c=='.'`; fold to lower (`dorp=='d'`) or upper. Digits, spaces, tabs, `@#~!=<>` are silently dropped. |
| `*` | `io.c:1380` `kake2hiku` | `*` → `-`. |
| Alphabet check | `mltaln9.c:60-85` `seqcheck`, called from `disttbfast.c:3470`, `tbfast.c:2483`, `pairlocalalign.c:3229`, `dvtditr.c:714`… | First char with `amino_n[c] == -1` → banner (via `reporterr`, silenced by `--quiet`) + `Illegal character c`, **exit 1**. Protein alphabet `blosum.c:12` `ARNDCQEGHILKMFPSTWYVBZX.-J`; nucleotide `DNA.h:43` `agctuAGCTUnNbdhkmnrsvwyx-O`. Hence protein **`U`, `O` are fatal**, `B Z X J` legal, **`.` is a legal protein residue** (index 23); nucleotide **`.` is fatal**, all IUPAC codes legal. |
| Gap removal | `mltaln9.c:10537` `gappick0` at `disttbfast.c:4450`, `tbfast.c:3320`, `pairlocalalign.c:3372`; `--add`: added `gappick0` `disttbfast.c:4315`, existing `commongappick` `disttbfast.c:4319` (`io.c:6296`) | Only **`-`** is removed. `.` is never removed. `--keeplength` does not change the stripping of added sequences. |

### 1b. `--anysymbol` / `--preservecase` (`anysymbol=1`, `mafft.tmpl:332-335`)

| Step | Where | Behaviour |
|---|---|---|
| Pre-pass | `mafft.tmpl:2305-2308` `replaceu $seqtype -i orig > infile` | Runs on the whole `infile`, which already includes the `--add` file (`mafft.tmpl:1142`). |
| Read | `replaceu.c:98,106` → `io.c:1634` `readData_pointer_casepreserve` → `io.c:1391` `load1SeqWithoutName_realloc_casepreserve` → `io.c:1329-1352` `charfilter` | Drops **only `\n`, space, `\r`**. Keeps digits, tabs, everything else. `=`, `<`, `>` inside a sequence → "Characters '= < >' can be used only in the title lines in the --anysymbol or --text mode." **exit 1**. |
| Substitute | `replaceu.c:129-137` | Protein: anything not in `ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.` → `X`, then `toupper`. Nucleotide: anything not in `ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-` → `n`, then `tolower`. Note `.` is usual for protein (stays a residue) but unusual for nucleotide (→ `n`). |
| Restore | `mafft.tmpl:2699` `restoreu -a pre -i orig`; `restoreu.c:170,203` `gappick_samestring` (`io.c:103-113`, removes `-` only), `restoreu.c:9-30` `fillorichar` | Original characters are copied onto every aligned position that is **not `-`**. Applies to added sequences too. |

## 2. Main vs fork B-2 vs C — behaviour comparison

| Sub-behaviour | C default | main (3e54a54) | fork B-2 |
|---|---|---|---|
| Drop digits/space/punct, fold case | yes (`onlyAlpha_*`) | yes (`normalize_sequence` + `apply_case_convention`) | yes + unconditional `replace_unusual` |
| `*` → `-` → removed | yes | yes | yes |
| Strip `-` from input before *every* phase | yes (`gappick0`) | **no on the L-INS pair path** (`engine.rs` pair-path `seq_refs` used the gapped input) | yes (`strip_input_gaps`) |
| `.` | protein: residue; DNA: `Illegal character .` exit 1 | stripped everywhere | **stripped** (invented; C never strips `.`) |
| `U`/`O` in protein, unknown letters | `Illegal character` exit 1 | aligned silently | **silently → `X`** (invented; C refuses) |
| `N` counts toward DNA detection | yes | **no** | no |
| Blank before `>` | exit 1 | accepted, header glued into previous sequence | accepted |
| `--anysymbol`: digits/tabs kept as residues | yes | **stripped** | stripped |
| `--anysymbol`: `=<>` in sequence | exit 1 | accepted | accepted |
| `--anysymbol` restore treats `.` | as residue | **as gap** (shifted restored residues) | as gap |
| `--anysymbol --add`: added seqs restored | yes | **no** (addfile read with the normalising reader) | no |

## 3. Corpus

`crates/mafft-bin/tests/fixtures/input_handling/` — 22 inputs (4 sequences of
~70 aa / ~65 bp each, C runs each in milliseconds), generated deterministically:

protein: `prot_lowercase`, `prot_star`, `prot_dot`, `prot_dash_inside`,
`prot_leading_trailing_gaps`, `prot_unusual_UOJBZX`, `prot_digits_spaces`
(GenBank-style numbering), `prot_crlf`, `prot_blank_lines`,
`prot_header_leading_ws`, `prot_punct` (`@ # ~ !` + tab), `prot_equals`;
DNA: `dna_mixed_case`, `dna_iupac` (N R Y K M S W B D H V), `dna_gaps_inside`
(`-` and `.`), `dna_leading_trailing_gaps`, `dna_crlf`, `dna_digits_spaces`,
`dna_n_rich` (every third base `N`), `dna_rna_u_star`;
`--add`: `add_existing_gapped.fa` (3 rows with a shared 6-column gap block) +
`add_new_unusual.fa` (lowercase + trailing `*`, and an internal `.`).

Modes: `default` (`--quiet`), `anysymbol` (`--anysymbol`), `linsi1`
(`--localpair --maxiterate 0`), `nuc` (DNA only, `--nuc`), and for `--add`:
`add`, `add_anysymbol`, `add_keeplength`. 71 (input, mode) cells. Every C
stdout is stored as `<name>.<mode>.expected`; exit-1 cases have no file and
are pinned by exit status + message in the tests.

## 4. Results

Legend: `=` byte-identical to C, `≠` differs, `E1` C exits 1 (stdout empty).
"main" = 3e54a54; "fixed" = this branch after the fix commit.

| input | mode | C | main | fixed | what differed on main |
|---|---|---|---|---|---|
| prot_lowercase | default/anysymbol/linsi1 | ok | = | = | |
| prot_star | default/anysymbol/linsi1 | ok | = | = | |
| prot_dot | default/anysymbol/linsi1 | ok | ≠ | = | `.` stripped instead of aligned as a residue |
| prot_dash_inside | default/anysymbol/linsi1 | ok | = | = | |
| prot_leading_trailing_gaps | default/anysymbol | ok | = | = | |
| prot_leading_trailing_gaps | linsi1 | ok | ≠ | = | pair phase saw the gapped input |
| prot_unusual_UOJBZX | default, linsi1 | E1 `Illegal character U` | exit 0 | E1 | main aligned `U`/`O` silently |
| prot_unusual_UOJBZX | anysymbol | ok | = | = | |
| prot_digits_spaces | default, linsi1 | ok | = | = | |
| prot_digits_spaces | anysymbol | ok | ≠ | = | C keeps the digits as `X` residues and restores them |
| prot_crlf | all | ok | = | = | |
| prot_blank_lines | all | ok | = | = | |
| prot_header_leading_ws | all | E1 (format message) | exit 0 | E1 | main glued `  >prot2` into `prot1` |
| prot_punct | default, linsi1 | ok | = | = | |
| prot_punct | anysymbol | ok | ≠ | = | tab kept by C as an `X` residue |
| prot_equals | default, linsi1 | ok | = | = | |
| prot_equals | anysymbol | E1 (`= < >` message) | exit 0 | E1 | |
| dna_mixed_case | all 4 | ok | = | = | |
| dna_iupac | all 4 | ok | = | = | |
| dna_gaps_inside | default, linsi1, nuc | E1 `Illegal character .` | exit 0 | E1 | |
| dna_gaps_inside | anysymbol | ok | ≠ | = | restore treated `.` as a gap → residues shifted |
| dna_leading_trailing_gaps | default/anysymbol/nuc | ok | = | = | |
| dna_leading_trailing_gaps | linsi1 | ok | ≠ | = | pair phase saw the gapped input |
| dna_crlf | all 4 | ok | = | = | |
| dna_digits_spaces | default, linsi1, nuc | ok | = | = | |
| dna_digits_spaces | anysymbol | ok | ≠ | = | digits kept by C as `n` residues |
| dna_n_rich | default, anysymbol, linsi1 | ok | ≠ | = | detected as protein (`N` not counted) |
| dna_n_rich | nuc | ok | = | = | |
| dna_rna_u_star | all 4 | ok | = | = | |
| add | add | ok | ≠ | = | `.` in the added sequence dropped |
| add | add_keeplength | ok | = | = | |
| add | add_anysymbol | ok | ≠ | **≠** | see residual R-A |

Summary: main diverged in 23 of 71 cells; fixed diverges in 1. Running the new
test file against main (fix files stashed): 13 of 24 tests fail, 11 pass,
1 ignored.

### Fork B-2 applied to main (same corpus)

Applying only the B-2 hunk (unconditional `replace_unusual` + `strip_input_gaps`)
to main and rerunning:

* fixes: `prot_leading_trailing_gaps.linsi1`, `dna_leading_trailing_gaps.linsi1`
  (the `-` strip reaches the pair phase);
* still wrong vs C: everything in the table involving `.`, `N`-detection,
  digits/tabs under `--anysymbol`, the two format errors, `--anysymbol --add`;
* newly wrong vs C: `prot_unusual_UOJBZX.default/linsi1` and
  `dna_gaps_inside.default/linsi1/nuc` now *succeed* with `U`→`X`, `.`→`n`
  where C exits 1 — and, because `strip_input_gaps` runs *after*
  `replace_unusual`, DNA `.` is turned into `n` and kept while protein `.`
  is deleted, which matches C in neither case.

## 5. Fixes made (separate commit, after the tests)

All are direct ports of the cited C behaviour; none is new policy.

1. `mafft-io/src/detect.rs` — `N` counts as nucleotide (`countATGC`).
2. `mafft-io/src/fasta.rs`, `error.rs` — case-preserving reader strips only `\n`,
   space, `\r` and rejects `= < >` (`charfilter`); both readers reject a
   blank-prefixed `>` line (`mafft.tmpl:1827`). Two new `IoError` variants carry
   C's messages.
3. `mafft-bin/src/lib.rs` — `seqcheck` after reading (and on the `--add` file):
   C's alphabets from `mafft_scoring`, C's banner unless `--quiet`, always
   `Illegal character c`, exit 1. Under `--anysymbol`, the `--add` file is read
   case-preserving, snapshotted and substituted like the main input, and the
   restore pass copies onto every non-`-` position (`fillorichar`).
   `read_error` passes C's two format-check messages through verbatim.
4. `mafft-core/src/engine.rs` — the pre-alignment strip removes `-` only, and the
   L-INS/G-INS/E-INS pair phase, `recompute_importance`, memsavetree and the
   RNA path all use that stripped set (`gappick0` in `pairlocalalign.c:3372`).
5. `mafft-core/src/add.rs`, `progressive.rs` — `.` is no longer treated as a gap
   in the `gappick0`/`commongappick` equivalents (C tests `== '-'` only). This
   is what made the `.` column vanish when two or more sequences were added.

### Residuals (not fixed, documented)

* **R-A** `add_anysymbol`: the added sequence ending in `*` (→ `X`) is placed
  `…KWRR-----X` by C and `…KWRRX-----` by rust. Reproduced with a literal
  trailing `X` and no `--anysymbol`, so it is a tie-break in the `--add`
  profile DP with a terminal zero-scoring residue, unrelated to input handling.
  Test `add_anysymbol_gapped_existing_with_unusual_new` is `#[ignore]`d with
  this reason.
* `mafft-align/src/profile.rs:207,336,480,503` still classify `.` as a gap when
  counting gap openings/closings. C tests `'-'` only
  (`mltaln9.c:12606-12660`). It changed no output in this corpus (protein `.`
  is rare), so it is left alone rather than touched blind.
* Under `--quiet`, C prints nothing at all for `Illegal character` (the message
  goes to `progressfile=/dev/null`); we keep the one-line message on stderr so
  a failing exit is never silent. stdout and exit status match.
* The blank-before-`>` message quotes the source line number; C quotes the line
  number in its CR-converted temp copy, which differs for CRLF input.

## 6. B-1 micro-benchmark

`crates/mafft-bin/examples/thread1_pool_bench.rs` — 8 worker threads, each
running `run_from(["--quiet","--thread","1","--retree","1","--maxiterate","0", fixture])`
N times; hashes every output against a reference. Wall time for
8 × 250 = 2000 calls, three runs each, same machine, back-to-back, release build:

| fixture | main (pool per call) | main + B-1 (cached thread-local pool) | saving |
|---|---|---|---|
| `input_handling/dna_mixed_case.fa` (4 × 65 bp) | 0.130–0.143 s | 0.106–0.120 s | ≈ 25 ms / 2000 calls ≈ **13 µs per call** of 8-way wall time (≈ 19 %) |
| `fixtures/refine_njob2_min.fa` (2 seqs) | 0.051–0.052 s | 0.040 s | ≈ 11 µs per call (≈ 22 %) |

Output identity: `mismatches=0` in every run, and `--thread 1` stdout for four
corpus files is `cmp`-identical between the two binaries. The saving is the
spawn/join of one rayon worker per call (~100 µs of thread time). An alignment
of realistic size takes 10²–10⁵ ms, so the gain is below noise there; it only
matters for a caller issuing thousands of sub-millisecond `--thread 1` calls
in-process.

Also noted while reading the hunk: the fork changes the per-call pool guard from
`args.thread > 0` to `> 1`, so `--thread 1` no longer builds a pool at all when
the cached one is unavailable — no output effect.

## 7. Recommendation

**B-1: accept, low priority.** Correct (thread-local, so no cross-thread
sharing; `OnceCell` so no double init), output-neutral, ~13 µs per call. Worth
taking as a small self-contained commit *if* pymafft/in-process users batch
many tiny alignments; otherwise cosmetic. Not applied on this branch — the
benchmark is committed so it can be re-measured after any pool refactor.

**B-2: reject as written.** Of its two parts:

* stripping `-` early is what C does (`gappick0`), and main was missing it on
  the pair path — fixed here at the engine level, where C does it, rather than
  in the CLI;
* stripping `.` and unconditionally applying `replace_unusual` are **invented**:
  C keeps protein `.` as a residue, refuses nucleotide `.`, and refuses
  protein `U`/`O` and any unknown letter with `Illegal character`. The fork
  would silently accept and rewrite inputs C rejects, which is exactly the
  class of divergence the byte-parity tests exist to catch.

What we must (and now do) match by default: `N` in type detection, `-`-only
gap stripping everywhere, `.`-as-protein-residue, `seqcheck` exit 1, the
blank-before-`>` exit 1; under `--anysymbol`: `charfilter`'s keep-everything
rule, the `= < >` exit 1, `-`-only restore, and restoring `--add` sequences.

## 8. Test runs (release, after the fix commit)

* `cargo test -p mafft-rs --release`: 50 + 5 + 24 passed, 0 failed, 1 ignored
  (the new `input_handling_vs_c.rs` contributes 24 passed / 1 ignored).
* `cargo test --workspace --exclude pymafft --release`: 493 passed, 0 failed, 7 ignored, exit 0.
* New test file against main (fix files stashed): 13 failed, 11 passed, 1 ignored.
