# MAFFT-rs

A Rust reimplementation of [MAFFT](https://mafft.cbrc.jp/alignment/software/), the widely-used multiple sequence alignment tool originally written in C by Kazutaka Katoh.

This project provides:
- **`mafft-rs`** — a single CLI binary replacing the original 5+ C binaries and shell script
- **`mafft-core`** — a Rust library crate for programmatic use
- **`pymafft`** — Python bindings via PyO3

## Status

**Production-ready.** On the included 36-sequence protein test dataset, the following modes produce **byte-identical output to C MAFFT 7.526** (`diff rust.fa c.fa` returns 0 lines):

- FFT-NS-2 (default), NW-NS-2 (`--nofft`)
- FFT-NS-i (`--maxiterate 100`)
- L-INS-1, G-INS-1, E-INS-1 (`--localpair` / `--globalpair` / `--genafpair` with `--maxiterate 0`)
- L-INS-i, G-INS-i, E-INS-i (with iterative refinement)
- BLOSUM 30 / 45 / 50 / 62 / 80 (FFT and NW)
- JTT 100 / 200, TM 100 / 200 (FFT and NW, including `--retree 1`)
- `--parttree`, `--dpparttree`, `--parttree --nofft` (divide-and-conquer guide tree)
- `--add`, `--add --nofft`, `--add --keeplength`
- Q-INS-i (RNA, requires `mxscarnamod`)
- `--allowshift`, `--reorder` / `--inputorder`, `--treeout`, `--treein`,
  `--auto`, `--anysymbol` / `--preservecase`,
  `--leavegappyregion` / `--legacygappenalty`, `--memsavetree`,
  `--memsave` / `--nomemsave` (Hirschberg `msalignmm` wired into engine),
  `--seed FILE`, `--seedtable FILE`, `--retree N` (clamped to [1,3]
  with INS-i forced to 1, matching `scripts/mafft:1840,1934`)

Every progressive merge step matches in score and width and every refinement iteration converges to C's exact alignment.

### Parity matrix (36-seq protein sample)

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
| BLOSUM 30 / 45 / 50 / 62 / 80 (FFT)         | match   | match      | 0    | ✓ byte-exact |
| BLOSUM 80 NW (`--bl 80 --nofft`)            | 712     | 712        | 0    | ✓ byte-exact |
| JTT 100 / 200 (FFT)                         | match   | match      | 0    | ✓ byte-exact |
| TM 100 / 200 (NW + FFT, default `--retree 2`) | match | match      | 0    | ✓ byte-exact |
| TM 200 (`--retree 1`)                       | 717     | 717        | 0    | ✓ byte-exact (closed 2026-06-02 — was previously misdiagnosed as §B.2 static-TLS; collateral from §B.14 FMA fusion fix) |
| `--parttree`, `--dpparttree`, `--parttree --nofft`, `--parttree --reorder` | 752 | 752 | 0 | ✓ byte-exact |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6 fixture) | match | match | 0 | ✓ byte-exact |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0‡   | ✓ byte-exact (case-insensitive) |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 0‡   | ✓ byte-exact (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0    | ✓ byte-exact |
| `--allowshift --globalpair --maxiterate 1000` | 1060   | 1060       | 0    | ✓ byte-exact (all 386 BALIBASE 3 files; X-strip ordering fix 2026-06-01) |
| `--reorder` (non-PartTree modes)            | match   | match      | 0    | ✓ byte-exact |

‡ For RNA / Q-INS-i: we uppercase residues; C preserves case. With `diff -i` RNA produces 0 lines.

Every supported mode is byte-identical to C MAFFT 7.526 on both the
36-seq protein test sample and the full BALIBASE 3 protein corpus
(1930/1930 alignments across 5 modes). Remaining items in `TODO.md`
are non-correctness gaps (performance deferrals and one
external-dependency limitation), not active divergences.

The original MAFFT C code is included as a git submodule for testing and cross-validation.

## Building

### Prerequisites

- Rust 1.80+ (edition 2024)
- Git (for submodule checkout)
- GCC/Clang (only needed to build the C reference for testing)

### Build

```bash
git clone --recurse-submodules https://github.com/luksgrin/rust-MAFFT.git
cd rust-MAFFT
cargo build --release
```

The binary is at `target/release/mafft-rs`.

### Run tests

```bash
# All integration tests across the workspace (357 tests, release build,
# pulls in mafft-sys for FFI cross-validation)
cargo test --workspace --release --tests

# Just unit tests (158 tests, faster; no mafft-sys / no C build)
cargo test --workspace --lib

# Just the end-to-end binary (89 tests, release build)
cargo test -p mafft-core --release --test end_to_end

# Build C reference for cross-validation (optional; cargo test --tests
# triggers this automatically via mafft-sys's build.rs)
make -C mafft-upstream/core
```

### Python bindings

```bash
cd crates/pymafft
uv venv .venv
uv pip install maturin pytest
maturin develop --release
uv run pytest tests/ -v   # 32 tests
```

## Usage

### CLI

```bash
# Basic usage (FFT-NS-2, fast default)
mafft-rs sequences.fasta > aligned.fasta

# Read from stdin
cat sequences.fasta | mafft-rs > aligned.fasta

# Write to file
mafft-rs -o aligned.fasta sequences.fasta

# Quiet mode (suppress progress messages)
mafft-rs -q sequences.fasta > aligned.fasta

# Control thread count (0 = all cores)
mafft-rs --thread 4 sequences.fasta > aligned.fasta
```

### Alignment strategies

```bash
# FFT-NS-2: fast progressive (default)
mafft-rs sequences.fasta

# FFT-NS-i: progressive + iterative refinement
mafft-rs --maxiterate 100 sequences.fasta

# G-INS-i: global iterative (accurate, for globally alignable sequences)
mafft-rs --globalpair --maxiterate 1000 sequences.fasta

# L-INS-i: local iterative (most accurate for <200 sequences)
mafft-rs --localpair --maxiterate 1000 sequences.fasta

# E-INS-i: generalized affine (for sequences with large internal gaps)
mafft-rs --genafpair --maxiterate 1000 sequences.fasta

# PartTree: divide-and-conquer guide tree (for 10K+ sequences)
mafft-rs --parttree sequences.fasta
```

### Adding sequences to an existing alignment

```bash
# Add new sequences to an existing aligned FASTA
mafft-rs existing_alignment.fasta --add new_sequences.fasta

# Add fragment sequences
mafft-rs existing_alignment.fasta --addfragments fragments.fasta

# Preserve existing alignment columns (no new gaps in original sequences)
mafft-rs existing_alignment.fasta --add new_sequences.fasta --keeplength
```

### Output formats

```bash
# FASTA (default)
mafft-rs sequences.fasta

# Clustal
mafft-rs --format clustal sequences.fasta

# PHYLIP
mafft-rs --format phylip sequences.fasta

# FASTA with custom line width (0 = no wrapping)
mafft-rs --linewidth 80 sequences.fasta
```

### Python

```python
import pymafft

# Align from list of sequences
result = pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])

# Align from named tuples
result = pymafft.align([("human", "ACDEFGHIK"), ("mouse", "ACDEFHIK")])

# Choose strategy
result = pymafft.align(seqs, strategy="linsi", maxiterate=1000)

# Align from file
result = pymafft.align_file("sequences.fasta")

# Access results
for seq in result:
    print(f"{seq.name}: {seq.sequence}")

result.to_fasta()    # FASTA string
result.to_tuples()   # list of (name, seq) tuples
```

## Migration from MAFFT (C)

The original MAFFT uses a shell script wrapper that invokes multiple C binaries. `mafft-rs` consolidates everything into a single binary with compatible flags.

### Command mapping

| Original MAFFT (C) | MAFFT-rs | Notes |
|---------------------|----------|-------|
| `mafft input.fa` | `mafft-rs input.fa` | FFT-NS-2 default |
| `mafft --maxiterate 100 input.fa` | `mafft-rs --maxiterate 100 input.fa` | FFT-NS-i |
| `mafft --globalpair input.fa` | `mafft-rs --globalpair input.fa` | G-INS-1 |
| `mafft --globalpair --maxiterate 1000 input.fa` | `mafft-rs --globalpair --maxiterate 1000 input.fa` | G-INS-i |
| `mafft --localpair input.fa` | `mafft-rs --localpair input.fa` | L-INS-1 |
| `mafft --localpair --maxiterate 1000 input.fa` | `mafft-rs --localpair --maxiterate 1000 input.fa` | L-INS-i |
| `mafft --ep 0 --genafpair --maxiterate 1000 input.fa` | `mafft-rs --genafpair --maxiterate 1000 input.fa` | E-INS-i |
| `linsi input.fa` | `mafft-rs --localpair --maxiterate 1000 input.fa` | L-INS-i shortcut |
| `ginsi input.fa` | `mafft-rs --globalpair --maxiterate 1000 input.fa` | G-INS-i shortcut |
| `einsi input.fa` | `mafft-rs --genafpair --maxiterate 1000 input.fa` | E-INS-i shortcut |
| `mafft --thread N input.fa` | `mafft-rs --thread N input.fa` | Multithreading |
| `mafft --op 1.53 input.fa` | `mafft-rs --op 1.53 input.fa` | Gap opening penalty |
| `mafft --ep 0.123 input.fa` | `mafft-rs --ep 0.123 input.fa` | Offset penalty |
| `mafft --bl 80 input.fa` | `mafft-rs --bl 80 input.fa` | BLOSUM matrix number |
| `mafft --jtt 200 input.fa` | `mafft-rs --jtt 200 input.fa` | JTT scoring matrix at PAM N |
| `mafft --tm 200 input.fa` | `mafft-rs --tm 200 input.fa` | Transmembrane scoring matrix at PAM N |
| `mafft --kimura 2 input.fa` | `mafft-rs --kimura 2 input.fa` | Kimura distance parameter |
| `mafft --allowshift --globalpair --maxiterate 1000 input.fa` | `mafft-rs --allowshift --globalpair --maxiterate 1000 input.fa` | G-INS-i with shift |
| `mafft --add new.fa existing.fa` | `mafft-rs existing.fa --add new.fa` | Add sequences to alignment |
| `mafft --addfragments frags.fa existing.fa` | `mafft-rs existing.fa --addfragments frags.fa` | Add fragment sequences |
| `mafft --nofft input.fa` | `mafft-rs --nofft input.fa` | Disable FFT (NW-NS-2) |
| `mafft --clustalout input.fa` | `mafft-rs --format clustal input.fa` | Clustal output |
| `mafft --phylipout input.fa` | `mafft-rs --format phylip input.fa` | PHYLIP output |

### Feature support status

| Original MAFFT flag | Status |
|---------------------|--------|
| `--maxiterate`, `--localpair`, `--globalpair`, `--genafpair` | Supported |
| `--add`, `--addfragments`, `--keeplength` | Supported |
| `--allowshift` | Supported (byte-identical to C) |
| `--nofft` | Supported |
| `--retree N` | Supported (clamp at 3, forces 1 for INS-i, matching `scripts/mafft:1840,1934`) |
| `--op`, `--ep`, `--bl` | Supported |
| `--thread N` | Supported (Rayon) |
| `--kimura N` | Supported (custom Kimura R for DNA PAM generation) |
| `--parttree`, `--dpparttree` | Supported (divide-and-conquer tree for 10K+ sequences) |
| `--groupsize N` | Supported (partition size for PartTree, default 150) |
| `--qinsi` (Q-INS-i) | Supported (requires `mxscarnamod` in PATH) |
| `--xinsi` (X-INS-i) | Supported (requires `contrafold` in PATH; untestable until installed) |
| `--scarnalike` | Supported (requires `dash_client` in PATH) |
| `--reorder` / `--inputorder` | Supported for all modes including PartTree |
| `--treeout` | Supported for all modes including `--parttree` / `--dpparttree` |
| `--memsave` / `--nomemsave` | Supported (Hirschberg `msalignmm` wired into engine for the non-FFT path) |
| `--seed`, `--seedtable`, `--treein`, `--auto`, `--anysymbol` / `--preservecase`, `--leavegappyregion` / `--legacygappenalty`, `--memsavetree` | Supported (byte-identical to C MAFFT 7.526) |

## Architecture

The project is organized as a Cargo workspace with 10 crates:

```
crates/
  mafft-sys/       Raw FFI bindings to MAFFT C code (dev-only, for cross-validation)
  mafft-types/     Shared Rust types (HomologyRegion, Sequence, ScoringContext, etc.)
  mafft-io/        FASTA, Clustal, PHYLIP, hat2 I/O
  mafft-scoring/   Substitution matrices (BLOSUM, JTT, TM, DNA) and gap penalties
  mafft-fft/       FFT-based homology detection (hand-ported Cooley-Tukey, bit-for-bit C-compat)
  mafft-align/     Pairwise and profile alignment algorithms (NW, SW, generalized affine, warp DP)
  mafft-tree/      Distance computation, NJ, UPGMA, guide tree construction, PartTree
  mafft-core/      Progressive alignment engine, iterative refinement, MafftEngine
  mafft-bin/       CLI binary (mafft-rs)
  pymafft/         Python bindings via PyO3
```

The release binary (`mafft-rs`) compiles with **zero C code** — `mafft-sys` is only used as a dev-dependency for cross-validation tests.

### Key design decisions

- **No global state.** The C code uses ~400 `extern` globals. Rust modules use owned `ScoringContext`, `Topology`, `Profile` structs passed explicitly.
- **`num_complex::Complex64`** replaces the C `Fukusosuu` struct. A hand-ported bit-for-bit Cooley-Tukey FFT in `mafft-fft/src/fft_c_compat.rs` matches C's `fft.c` rounding exactly (replaced `rustfft` after BL50 / TM 200 FFT divergences were traced to 1-ULP correlation differences).
- **`f64` DP matrices** (post 2026-05-10 migration). The substitution matrix is stored as `Vec<Vec<i32>>` for canonical storage and exposed as `consweight_matrix: Vec<Vec<f64>>` for DP arithmetic, preserving sub-integer precision through per-cell scoring (matches C's `n_dis_consweight_multi`).
- **`f64::mul_add` (FMA) in DP arithmetic** — gcc `-O3` fuses `a + b * c` into a single-rounding FMA, while plain Rust `+=` produces two roundings. Rust uses `mul_add` in `match_calc_row` and gap-frequency-modulated penalties to match C bit-for-bit on flat-landscape matrices (closed the BL50 divergence).
- **Rayon parallelism** for pairwise distance computation, all-vs-all local alignments, and refinement scoring. Sequential float summation is preserved for accept/reject decisions to keep refinement output deterministic.
- **SIMD-friendly inner loops**: branchless patterns in `match_score()`, `pairwise_score()`, and `pairwise_identity_distance()` that LLVM auto-vectorizes to NEON/AVX instructions.
- **Three-way gap insertion** (matching C's `insertnewgaps()`): group1 follows cursor1, group2 follows cursor2, "other" sequences follow cursor1 with gaps at Insert positions.
- **Platform-independent sort tie-break.** C MAFFT calls `qsort()` directly inside `splitseq_mq` to sort sequences by their distance to the pivot. `qsort` is part of the host C library — BSD on macOS, glibc on Linux, MSVC on Windows — and its three implementations disagree on the relative order of *truly tied* elements (same distance, same selfscore, same length, which only happens when the input contains exactly-duplicate sequences). That makes C MAFFT's `--parttree --reorder` output platform-dependent: the same C source compiled on macOS vs Linux produces different orderings for duplicates. Our Rust port ships a hand-written BSD-qsort algorithm (`mafft-tree/src/bsd_qsort.rs`, Bentley-McIlroy) and uses it everywhere `dcompare_sort` is called, so the `mafft-rs` binary produces the **same output on every platform it's built for**, matching macOS C MAFFT 7.526 byte-for-byte.

### Test suite

`cargo test --workspace --release --tests`: **364 passed, 0 failed,
4 ignored** (4 are deliberate diagnostic/bisection tools, annotated
at the call site). Lib tests separately: `cargo test --workspace
--lib`: **158 passed**. Python: 32 tests pass.

| Suite | Count | What |
|-------|-------|------|
| Rust lib tests (`--lib`, all crates) | 158 | All `#[cfg(test)] mod tests` blocks (incl. 7 `msalign` Hirschberg DP tests + 4 `parse_hat3_seed` tests) |
| Rust binary tests (`mafft-rs`) | 10 | CLI helper functions (`decide_auto`, `replace_unusual`, etc.) |
| Rust integration tests (`end_to_end`) | 89 | Byte-level parity with C across all supported modes + DP diagnostics (incl. 4 `--seedtable` and 7 `--retree N` byte-equal tests + the BB30013 boundary-frequencies regression guard) |
| Rust FFI cross-validation tests | ~190 | Cell-by-cell matrix equality, single-pair `G__align11` / `A__align` / `genL__align11` / warp DP / FFT equality, PartTree pipeline, memsavetree, `MSalignmm` Hirschberg DP, MSalignmm profile alignment, BB20027 weight-equality, BB12041 distance-equality |
| Rust other integration tests | ~25 | `trace_refinement` (15), `integration` mafft-io (7), `imp_*` (3) |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total Rust integration** | **364** | |

Regression guards for C parity are in `crates/mafft-core/tests/end_to_end.rs` (mode-level byte-identity), `crates/mafft-core/tests/cross_validate_*.rs` (FFI-level cell/function equality), and `crates/mafft-tree/tests/cross_validate_parttree.rs` (PartTree pipeline equality). Any regression in DP indexing, boundary handling, FFT anchor segment gaps, retree distance, refinement-tree distance, or pairwise/profile consistency will fail at least one of these.

#### BALIBASE 3 parity sweep

Sweep across BALIBASE 3 (RV11–RV50, 386 protein test sets) against MAFFT 7.526:

| Mode | Flags | Match | Total % | Residuals |
|------|-------|-------|---------|-----------|
| FFT-NS-i | `--maxiterate 1000` | **386/386** | **100.0 %** | none |
| L-INS-i  | `--localpair --maxiterate 1000` | **386/386** | **100.0 %** | none |
| G-INS-i  | `--globalpair --maxiterate 1000` | **386/386** | **100.0 %** | none |
| E-INS-i  | `--genafpair --maxiterate 1000` | **386/386** | **100.0 %** | none |
| **G-INS-i + allowshift** | `--globalpair --maxiterate 1000 --allowshift` | **386/386** | **100.0 %** | none (2026-06-01) |

**All five sweep modes — including the previously-problematic
`--allowshift` G-INS-i — are now 386/386 byte-identical.** Total:
**1930/1930 (100 %)** across the BALIBASE 3 protein corpus.

The `--allowshift` residuals (originally 45 of 386, then 18 after
performance work) closed via the X-strip ordering fix
(`constraints.rs::build_homology_table_with_unalign`,
2026-06-01). C's `pairlocalalign.c:2197` runs
`G__align11_noalign` on X-stripped sequences *before* computing the
distance that feeds `makedynamicmtx` (line 2207–2215). Rust was doing
the X-strip recompute *after* the dynmtx re-run, so for X-containing
inputs the dynamic matrix was built from the wrong distance →
different second-alignment regions → different LH `importance` →
final-alignment flips. Fix: compute `pscore_for_dist` upfront,
apply `alignment.score = pscore_for_dist` unconditionally for L/G
aligners. This also definitively disproved the previously-suspected
"C `A__align` static-TLS state machine" hypothesis: an instrumented C
binary that **zeros every static-TLS buffer at A__align entry**
produces byte-identical output to vanilla C. The 45 §B.2-class
"tied-trace" cases were a single FP-input-ordering bug all along, not
unfixable C-side state.

The earlier four INS-i residuals (BB30028, BB50001 L-INS-i; BB20004
G-INS-i; BB40004 E-INS-i) closed via FP summation-order and tree-phase
fixes in the constrained refinement path. None was the static-TLS
artifact once suspected. All were found by instrumenting C MAFFT and
comparing bit patterns at `%.18e`:

1. **`hat2` 3-decimal truncation** for the refinement-weights tree
   (`engine.rs`) — C's `dvtditr` rebuilds the tree from the 3-decimal
   `hat2` file, not the in-memory full-precision matrix. Closed
   BB50001 (244-line residual).
2. **printf banker's rounding for `opt`** (`constraints.rs`) — C's
   `%7.5f` hat3 serialization rounds half-to-even; Rust's `f64::round`
   rounds half-away-from-zero. At exact half-points (`x.x005`) they
   disagree by 1 ULP. `format!("{:.5}", v).parse()` matches C, making
   the importance table bit-identical (12694/12694 regions).
3. **backward diagonal sum for `old_imp`** (`refinement.rs`) — C's
   `tditeration.c:891` reads `imp_match_out_scD` from `i=length-1`
   down; FP addition is non-associative, so direction matters. Closed
   BB20004 (10-line residual).
4. **per-segment forward `impmatch` accumulation for `new_imp`**
   (`refinement.rs`) — C's `Falign_localhom` sums each FFT segment's
   impmatch (backward within a segment) then accumulates forward
   across segments (`Falign_localhom.c:816`); a single global diagonal
   sweep diverges by ~7 ULP. Closed BB30028 (4-line residual).
5. **phase-split constraint importance** (`engine.rs`) — C computes the
   constraint `importance` twice with *different* trees: tbfast
   (progressive) uses the full-precision in-memory `iscore` tree, while
   dvtditr (refinement) re-reads the 3-decimal `hat2` tree. Rust used
   the rounded tree for both. The progressive phase now uses the
   full-precision tree; a second `recompute_importance` from the rounded
   refinement tree runs before `iterative_refine`. The E-INS-i BB40004
   divergence was purely this (the genaffine alignment + distances were
   already bit-identical; `--treein` with C's tree gave diff=0). Closed
   BB40004 (2901-line residual); the per-phase split also makes L/G-INS-i
   robustly correct rather than coincidentally matching.

Cell-level FFI tests (`cross_validate_cpmx`, `cross_validate_counteff`,
`cross_validate_bb20027_dp`) prove the Rust DP and profile-building
match C bit-for-bit on identical inputs.

FFT-NS-i reached 100 % earlier via the §B.10 FFT-segmented refinement
boundary-frequencies fix (2026-05-25): closed BB30013's 3843-line
residual by computing per-segment `headgapfreq{1,2}` /
`gapfreq{1,2}[lgth]` from `sgap`/`egap` via `outgapcount` (mirroring
C `Salignmm.c:1585-1622`) instead of defaulting `BoundaryFreqs` to
1.0. Full writeup in `TODO.md §B.10`.

To reproduce:

```bash
cargo build --release
scripts/balibase_parity.py /path/to/BAliBASE/ \
    --modes "" "--maxiterate 100" "--localpair --maxiterate 100" \
    --output balibase_parity.tsv
```

Full investigation log: `balibase_parity_run.md`.

## Known limitations

### `--xinsi` (untestable)

Wired correctly but requires Stanford's `CONTRAfold v2.02+` binary, which is not shipped by upstream MAFFT. The binary emits a "contrafold not found" diagnostic when run without the external dependency. See `TODO.md` §D.

### Missing CLI flags

None. Every user-visible MAFFT 7.526 CLI flag is implemented and
byte-identical to C on the 36-seq test sample AND across all 386
BALIBASE 3 protein test sets (verified for FFT-NS-i, L-INS-i, G-INS-i,
E-INS-i, and `--allowshift` G-INS-i, all at `--maxiterate 1000`).

The remaining `TODO.md` items are:
- `--xinsi` untestable until Stanford's `contrafold` is installed
  (§D — wired and emits the correct "contrafold not found" diagnostic).

**Zero known correctness divergences** across the 36-seq protein test
sample and the BALIBASE 3 protein corpus (1930/1930 alignments matched
byte-for-byte to C MAFFT 7.526).

### Case preservation

C preserves the input case of residues (e.g., lowercase RNA). We uppercase all residues before alignment. The alignment itself (gap placement) is identical — `rna_nofft_case_insensitive_identical_to_c` verifies this.

### Performance

Benchmark snapshot vs C MAFFT 7.526 (2026-06-02, macOS arm64):

| Mode / input                                | C MAFFT | mafft-rs | Ratio              |
|---------------------------------------------|---------|----------|--------------------|
| default FFT-NS-2 (108-seq synthetic)        | 1.34s   | 1.16s    | **rust 1.16× faster** |
| `--maxiterate 100` FFT-NS-i (36-seq sample) | 0.767s  | 0.737s   | **rust 1.04× faster** |
| `--maxiterate 100` FFT-NS-i (BB12019)       | 0.302s  | 0.205s   | **rust 1.47× faster** |
| `--maxiterate 100` FFT-NS-i (BB30004)       | 1.370s  | 1.384s   | rust 0.99× (within noise) |
| `--allowshift --globalpair --maxiterate 1000` (36-seq) | 2.56s | 1.73s | **rust 1.48× faster** |
| `--maxiterate 50 --localpair` (108-seq)     | 37.00s  | 22.93s   | **rust 1.61× faster** |
| `--maxiterate 50 --globalpair` (108-seq)    | 35.32s  | 22.62s   | **rust 1.56× faster** |
| `--parttree` (108-seq)                      | 1.70s   | 0.09s    | **rust 18× faster**   |
| `--dpparttree` (108-seq)                    | 7.24s   | 0.08s    | **rust 90× faster**   |

FFT-NS-i closed 2026-06-02 — Rust is now at-or-better than C across the
inputs we measure. The closing change was a focused pass over
`profile_align_imp_multimtx`: `unsafe { get_unchecked }` on every hot
inner-loop read/write, hoisting of two loop-invariants (`fgcp1[i-1]`,
`ogcp1[i]`) out of the j-loop, slice-binding of the per-row-constant
arrays to fixed-lifetime locals, and gating the six warp-state buffer
allocations behind `try_warp` (so FFT-NS-i no longer pays for warp
setup it never uses; the same gate turned out to be a 1.48× win on
`--allowshift` itself by trimming setup time of the warp DP path too).
Each `get_unchecked` call site is annotated with the SAFETY invariant
that makes it sound. Full details in `TODO.md §C.2`; byte-identity to
C MAFFT 7.526 was preserved at every step (90/90 end-to-end +
15/15 cross-validate FFI cell-equality tests pass; direct binary diff
is 0 on the 36-seq sample for both FFT-NS-i and `--allowshift
--globalpair --maxiterate 1000`, plus BB30004 and BB12019).

Earlier-flagged §C.1 (per-group gap stripping) is documented as a
**design choice** rather than deferred work — the modes it would have
targeted (FFT-NS-2, L/G-INS-i, PartTree, DP-PartTree) are already
1.16×–90× faster than C without it.

SIMD-friendly patterns in `match_score()`, `pairwise_score()`,
`pairwise_identity_distance()` auto-vectorize via LLVM.

## Upstream MAFFT

The original MAFFT C code is included as a git submodule at `mafft-upstream/`, pinned to commit `ee97999` (version 7.526). To update:

```bash
cd mafft-upstream
git fetch origin
git checkout <new-commit-or-tag>
cd ..
git add mafft-upstream
cargo test --workspace  # verify nothing breaks
git commit -m "Bump MAFFT upstream to <version>"
```

## CI/CD

| Workflow | Triggers | What |
|----------|----------|------|
| **CI** (`ci.yml`) | push/PR to main/dev | Build + test Rust and C, smoke test binary |
| **Python** (`python.yml`) | push/PR + release | Build wheels (Linux, macOS, Windows), test on Python 3.9-3.13, publish to PyPI on release |
| **Release** (`release.yml`) | GitHub release | Build `mafft-rs` binaries for 4 platforms, attach to release |

## License

- Rust code (`crates/`): [MIT](LICENSE-MIT)
- Original MAFFT C code (`mafft-upstream/`): [BSD-3-Clause](mafft-upstream/license)

## Credits

- Original MAFFT by [Kazutaka Katoh](https://mafft.cbrc.jp/alignment/software/) (CBRC, AIST)
- Rust reimplementation by Lucas Goiriz (luksgrin)
