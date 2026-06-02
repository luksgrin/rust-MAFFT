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
| TM 200 (`--retree 1`)                       | 717     | 717        | 0    | ✓ byte-exact |
| `--parttree`, `--dpparttree`, `--parttree --nofft`, `--parttree --reorder` | 752 | 752 | 0 | ✓ byte-exact |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6 fixture) | match | match | 0 | ✓ byte-exact |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0‡   | ✓ byte-exact (case-insensitive) |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 0‡   | ✓ byte-exact (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0    | ✓ byte-exact |
| `--allowshift --globalpair --maxiterate 1000` | 1060   | 1060       | 0    | ✓ byte-exact (all 386 BALIBASE 3 files) |
| `--reorder` (non-PartTree modes)            | match   | match      | 0    | ✓ byte-exact |

‡ For RNA / Q-INS-i: we uppercase residues; C preserves case. With `diff -i` RNA produces 0 lines.

Every supported mode is byte-identical to C MAFFT 7.526 on both the
36-seq protein test sample and the full BALIBASE 3 protein corpus
(1930/1930 alignments across 5 modes). Remaining items in `TODO.md`
are non-correctness gaps (unwired CLI flags, external-dep-blocked
modes), not active divergences.

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

The binary is at `target/release/mafft-rs`. To match C MAFFT's
symlink-based shortcuts (`linsi` / `ginsi` / `einsi` / `fftns` /
`fftnsi` / `nwns` / `nwnsi` / `qinsi` / `xinsi`), drop the binary on
your `PATH` and create symlinks:

```bash
ln -s mafft-rs linsi
ln -s mafft-rs ginsi
ln -s mafft-rs einsi
# ... etc.
```

`mafft-rs` reads its `argv[0]` basename (stripping an optional
`mafft-` prefix) and applies the same defaults as C's shell script.
Explicit flags still override (e.g. `linsi --maxiterate 5 ...`
honours the override).

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
| `linsi input.fa` | `linsi input.fa` (symlink) or `mafft-rs --localpair --maxiterate 1000 input.fa` | L-INS-i shortcut |
| `ginsi input.fa` | `ginsi input.fa` (symlink) or `mafft-rs --globalpair --maxiterate 1000 input.fa` | G-INS-i shortcut |
| `einsi input.fa` | `einsi input.fa` (symlink) or `mafft-rs --genafpair --maxiterate 1000 input.fa` | E-INS-i shortcut |
| `fftnsi input.fa` | `fftnsi input.fa` (symlink) or `mafft-rs --maxiterate 2 input.fa` | FFT-NS-i shortcut (C default is `--maxiterate 2`, not 100) |
| `nwns input.fa` | `nwns input.fa` (symlink) or `mafft-rs --nofft input.fa` | NW-NS-2 shortcut |
| `nwnsi input.fa` | `nwnsi input.fa` (symlink) or `mafft-rs --nofft --maxiterate 2 input.fa` | NW-NS-i shortcut |
| `qinsi input.fa` | `qinsi input.fa` (symlink) or `mafft-rs --qinsi --maxiterate 1000 input.fa` | Q-INS-i shortcut |
| `xinsi input.fa` | `xinsi input.fa` (symlink) or `mafft-rs --xinsi --maxiterate 1000 input.fa` | X-INS-i shortcut |
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
- **`num_complex::Complex64`** replaces the C `Fukusosuu` struct. A hand-ported bit-for-bit Cooley-Tukey FFT in `mafft-fft/src/fft_c_compat.rs` matches C's `fft.c` rounding exactly. (Necessary because off-the-shelf FFT libraries vary in butterfly grouping order, which produces 1-ULP correlation differences that flip FFT-anchor selection on flat-landscape matrices.)
- **`f64` DP matrices.** The substitution matrix is stored as `Vec<Vec<i32>>` for canonical storage and exposed as `consweight_matrix: Vec<Vec<f64>>` for DP arithmetic, preserving sub-integer precision through per-cell scoring (matches C's `n_dis_consweight_multi`).
- **`f64::mul_add` (FMA) in DP arithmetic.** gcc / clang `-O3` with `FP_CONTRACT=on` fuse `a + b * c` into a single-rounding FMA; plain Rust `+=` produces two roundings. Rust uses `mul_add` in `match_calc_row` and gap-frequency-modulated penalties to match C bit-for-bit on flat-landscape matrices.
- **Rayon parallelism** for pairwise distance computation, all-vs-all local alignments, and refinement scoring. Sequential float summation is preserved for accept/reject decisions to keep refinement output deterministic.
- **SIMD-friendly inner loops**: branchless patterns in `match_score()`, `pairwise_score()`, and `pairwise_identity_distance()` that LLVM auto-vectorizes to NEON/AVX instructions.
- **Three-way gap insertion** (matching C's `insertnewgaps()`): group1 follows cursor1, group2 follows cursor2, "other" sequences follow cursor1 with gaps at Insert positions.
- **Platform-independent sort tie-break.** C MAFFT calls `qsort()` directly inside `splitseq_mq` to sort sequences by their distance to the pivot. `qsort` is part of the host C library — BSD on macOS, glibc on Linux, MSVC on Windows — and its three implementations disagree on the relative order of *truly tied* elements (same distance, same selfscore, same length, which only happens when the input contains exactly-duplicate sequences). That makes C MAFFT's `--parttree --reorder` output platform-dependent: the same C source compiled on macOS vs Linux produces different orderings for duplicates. Our Rust port ships a hand-written BSD-qsort algorithm (`mafft-tree/src/bsd_qsort.rs`, Bentley-McIlroy) and uses it everywhere `dcompare_sort` is called, so the `mafft-rs` binary produces the **same output on every platform it's built for**, matching macOS C MAFFT 7.526 byte-for-byte.

### Test suite

`cargo test --workspace --release --tests -- --test-threads=1`:
**364 passed, 0 failed, 4 ignored**. The 4 ignored tests are
deliberate diagnostic / bisection tools annotated at the call site.
Python (`pymafft`): **32 passed**.

| Suite | Count | What |
|-------|-------|------|
| Rust lib tests (`--lib`, all crates)  | 158  | All `#[cfg(test)] mod tests` blocks (incl. 7 `msalign` Hirschberg DP tests + 4 `parse_hat3_seed` tests) |
| Rust binary tests (`mafft-rs`)        | 10   | CLI helper functions (`decide_auto`, `replace_unusual`, etc.) |
| Rust integration (`end_to_end`)       | 90   | Mode-level byte-identity vs C across every supported flag/mode combo |
| Rust FFI cross-validation             | ~90  | Cell-by-cell FFI equality: `G__align11`, `A__align`, `genL__align11`, warp DP, FFT, MSalignmm, PartTree pipeline, memsavetree, BB20027 weights, BB12041 distances |
| Other integration                     | ~16  | `trace_refinement`, `integration` (mafft-io), `imp_*` |
| Python tests                          | 32   | API, strategies, file I/O, error handling, types |
| **Total Rust**                        | **364** | |

Regression guards for C parity live in
`crates/mafft-core/tests/end_to_end.rs` (mode-level byte-identity),
`crates/mafft-core/tests/cross_validate_*.rs` (FFI-level cell/function
equality), and `crates/mafft-tree/tests/cross_validate_parttree.rs`
(PartTree pipeline equality). Any regression in DP indexing, boundary
handling, FFT anchor segment gaps, retree distance,
refinement-tree distance, or pairwise/profile consistency fails at
least one of them.

**Serialise the workspace test run** (`--test-threads=1`). One of the
FFI cross-validation suites has a teardown-time `free()` issue under
parallel execution; serialising avoids it. Individual suites run fine
parallel; only the full-workspace combine triggers the flake.

#### BALIBASE 3 parity sweep

Sweep across BALIBASE 3 (RV11–RV50, 386 protein test sets) against MAFFT 7.526:

| Mode | Flags | Match | Total % | Residuals |
|------|-------|-------|---------|-----------|
| FFT-NS-i | `--maxiterate 1000` | **386/386** | **100.0 %** | none |
| L-INS-i  | `--localpair --maxiterate 1000` | **386/386** | **100.0 %** | none |
| G-INS-i  | `--globalpair --maxiterate 1000` | **386/386** | **100.0 %** | none |
| E-INS-i  | `--genafpair --maxiterate 1000` | **386/386** | **100.0 %** | none |
| **G-INS-i + allowshift** | `--globalpair --maxiterate 1000 --allowshift` | **386/386** | **100.0 %** | none (2026-06-01) |

Total: **1930/1930 alignments byte-identical** to C MAFFT 7.526
across the BALIBASE 3 protein corpus. Cell-level FFI tests
(`cross_validate_cpmx`, `cross_validate_counteff`,
`cross_validate_bb20027_dp`) additionally prove the Rust DP and
profile-building match C bit-for-bit on identical inputs.

To reproduce:

```bash
cargo build --release
scripts/balibase_parity.py /path/to/BAliBASE/ \
    --modes "" "--maxiterate 100" "--localpair --maxiterate 100" \
    --output balibase_parity.tsv
```

## Known limitations

### External-dependency-blocked modes

These three modes are wired correctly but require third-party binaries
not shipped by upstream MAFFT; end-to-end validation depends on the
binary being on `PATH`:

- **`--xinsi` (X-INS-i)** — needs Stanford's
  [CONTRAfold v2.02+](http://contra.stanford.edu/contrafold/). Emits a
  "contrafold not found" diagnostic if absent.
- **`--qinsi` (Q-INS-i)** — needs `mxscarnamod` built from
  `mafft-upstream/extensions`.
- **`--scarnalike`** — needs `dash_client` in `PATH`.

### Unsupported CLI flags

The Rust binary covers every CLI flag exercised by the BALIBASE 3
sweep and every flag in the parity matrix above. A number of niche C
MAFFT 7.526 flags are NOT implemented; see `TODO.md` "Active items"
for the full inventory. Notable omissions:

- **No `argv[0]` shortcut dispatch.** C MAFFT ships `linsi` / `ginsi`
  / `einsi` / `fftnsi` / `nwnsi` / `qinsi` / `xinsi` symlinks; `mafft-
  rs` requires long-flag form. See "Migration from MAFFT (C)" below
  for the equivalent `mafft-rs` invocations.
- **No DNA strand auto-detection** (`--adjustdirection`,
  `--adjustdirectionaccurately`).
- **No alternative tree-linkage flags** (`--averagelinkage`,
  `--minimumlinkage`, `--mixedlinkage`, `--youngestlinkage`); the C
  default is used.
- **No alternative iteration-strategy flags** (`--bestfirst`,
  `--simplehillclimbing`, `--skipiterate`, `--oneiteration`).
- **No auxiliary output flags** beyond `--format clustal/phylip` and
  `--treeout` (`--distout`, `--scoreout`, `--mapout`, `--pileup`,
  etc. are not implemented).
- **No RNA-structure or PDB-structure alternative pipelines**
  (DAFS, FoldAlign, LARA, SCARNA, MCCASKILL, RIBOSUM, RNAALIFOLD,
  `--pdbfilelist`, `--pdbidlist`) — these are out of scope.
- **No MPI flag** (`--mpi`) — we use Rayon.

### Case preservation

C preserves the input case of residues (e.g., lowercase RNA). We uppercase all residues before alignment. The alignment itself (gap placement) is identical — `rna_nofft_case_insensitive_identical_to_c` verifies this.

### Performance

Benchmark snapshot vs C MAFFT 7.526 (macOS arm64, alternating-run medians):

| Mode / input                                | C MAFFT | mafft-rs | Ratio              |
|---------------------------------------------|---------|----------|--------------------|
| default FFT-NS-2 (108-seq synthetic)        | 1.34s   | 1.16s    | **rust 1.16× faster** |
| `--maxiterate 100` FFT-NS-i (36-seq sample) | 0.785s  | 0.670s   | **rust 1.17× faster** |
| `--maxiterate 100` FFT-NS-i (BB12019)       | 0.333s  | 0.198s   | **rust 1.68× faster** |
| `--maxiterate 100` FFT-NS-i (BB30004)       | 1.366s  | 1.246s   | **rust 1.10× faster** |
| `--allowshift --globalpair --maxiterate 1000` (36-seq) | 2.586s | 1.756s | **rust 1.47× faster** |
| `--maxiterate 50 --localpair` (108-seq)     | 37.00s  | 22.93s   | **rust 1.61× faster** |
| `--maxiterate 50 --globalpair` (108-seq)    | 35.32s  | 22.62s   | **rust 1.56× faster** |
| `--parttree` (108-seq)                      | 1.70s   | 0.09s    | **rust 18× faster**   |
| `--dpparttree` (108-seq)                    | 7.24s   | 0.08s    | **rust 90× faster**   |

SIMD-friendly patterns in `match_score()`, `pairwise_score()`, and
`pairwise_identity_distance()` auto-vectorize via LLVM. The hot inner
DP loop (`profile_align_imp_multimtx` in
`crates/mafft-align/src/profile.rs`) is specialised at codegen time
via a `#[inline(always)] fn j_loop<const TW: bool, const STRICT:
bool>` helper, dispatched from the i-loop with a 4-way match — LLVM
strips the warp arithmetic from the FFT-NS-i path and replaces
per-cell branches on `strict_part_tiebreak` with const choices. All
per-call scratch vectors come from thread-local pools (`DpScratch` /
`DP_H_POOL` / `DP_IJP_POOL`) to amortise allocation across calls,
matching C `A__align`'s static-TLS buffer strategy.

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
