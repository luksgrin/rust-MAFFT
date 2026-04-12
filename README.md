# MAFFT-rs

A Rust reimplementation of [MAFFT](https://mafft.cbrc.jp/alignment/software/), the widely-used multiple sequence alignment tool originally written in C by Kazutaka Katoh.

This project provides:
- **`mafft-rs`** — a single CLI binary replacing the original 5+ C binaries and shell script
- **`mafft-core`** — a Rust library crate for programmatic use
- **`pymafft`** — Python bindings via PyO3

## Status

**Working implementation.** The core alignment pipeline (progressive alignment, iterative refinement, FFT-accelerated homology detection) is implemented and produces valid alignments for protein and DNA sequences. On the included 36-sequence test dataset, alignment quality exceeds the C implementation (131% SP ratio). See [Gap Analysis](#gap-analysis-vs-original-mafft) for remaining work.

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
# Fast unit tests (all crates, ~91 tests)
cargo test --workspace --lib

# Integration tests (requires release build, ~4 tests)
cargo test -p mafft-core --release --test end_to_end

# Build and test C reference (optional)
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
| `--allowshift` | Supported |
| `--nofft` | Supported |
| `--retree N` | Supported (default 2, matching C) |
| `--op`, `--ep`, `--bl` | Supported |
| `--thread N` | Supported (Rayon) |
| `--kimura N` | Supported (custom Kimura R for DNA PAM generation) |
| `--parttree`, `--dpparttree` | Supported (divide-and-conquer tree for 10K+ sequences) |
| `--groupsize N` | Supported (partition size for PartTree, default 150) |
| `--qinsi` (Q-INS-i) | Supported (requires `mxscarnamod` in PATH) |
| `--xinsi` (X-INS-i) | Supported (requires `contrafold` in PATH) |
| `--scarnalike` | Supported (requires `dash_client` in PATH) |

## Architecture

The project is organized as a Cargo workspace with 10 crates:

```
crates/
  mafft-sys/       Raw FFI bindings to MAFFT C code (dev-only, for cross-validation)
  mafft-types/     Shared Rust types (HomologyRegion, Sequence, ScoringContext, etc.)
  mafft-io/        FASTA, Clustal, PHYLIP, hat2 I/O
  mafft-scoring/   Substitution matrices (BLOSUM, JTT, TM, DNA) and gap penalties
  mafft-fft/       FFT-based homology detection (using rustfft + num-complex)
  mafft-align/     Pairwise and profile alignment algorithms (NW, SW, generalized affine)
  mafft-tree/      Distance computation, NJ, UPGMA, guide tree construction
  mafft-core/      Progressive alignment engine, iterative refinement, MafftEngine
  mafft-bin/       CLI binary (mafft-rs)
  pymafft/         Python bindings via PyO3
```

The release binary (`mafft-rs`) compiles with **zero C code** — `mafft-sys` is only used as a dev-dependency for cross-validation tests.

### Key design decisions

- **No global state.** The C code uses ~400 `extern` globals. Rust modules use owned `ScoringContext`, `Topology`, `Profile` structs passed explicitly.
- **`num_complex::Complex64`** replaces the C `Fukusosuu` struct. `rustfft` replaces the hand-rolled Cooley-Tukey FFT.
- **Rayon parallelism** for pairwise distance computation, all-vs-all local alignments, and refinement scoring.
- **SIMD-friendly inner loops**: branchless patterns in `match_score()`, `pairwise_score()`, and `pairwise_identity_distance()` that LLVM auto-vectorizes to NEON/AVX instructions.
- **Three-way gap insertion** (matching C's `insertnewgaps()`): group1 follows cursor1, group2 follows cursor2, "other" sequences follow cursor1 with gaps at Insert positions.

### Test suite

| Suite | Count | What |
|-------|-------|------|
| Rust unit tests | 103 | All crates, all modules |
| Rust integration tests | 4 | End-to-end on real data + C reference comparison |
| C alignment tests | 8 | FFT-NS-2, FFT-NS-i, G-INS-i, L-INS-i, parttree, etc. |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total** | **148** | |

## Known limitations

### Alignment quality

On the included 36-sequence protein test dataset (`mafft-upstream/test/sample`), the Rust implementation achieves **97.0%** of the C implementation's sum-of-pairs identity score (SP=0.316 vs C's 0.326, width 703 vs 717) for FFT-NS-2. The profile alignment DP now ports C's exact formulation: gap opening from diagonal, substitution score on all paths, `ijp` traceback with jump distances. The alignment is structurally correct (all sequences have the same width, ungapped sequences match originals, residue content is preserved).

### Performance

- No SIMD for the DP fill loops themselves (data dependencies prevent vectorization without anti-diagonal restructuring).

## Remaining gaps vs original MAFFT

### Alignment quality (3.3% SP score gap)

The profile alignment DP ports C's exact formulation (diagonal gap opening, sub on all paths, `ijp` traceback, `st_OpeningGapCount`/`st_FinalGapCount`, headgapfreq=1.0, profile caching). The remaining gap is from:

- **Profile caching granularity**: C's `cpmxhist` blends profiles INSIDE the rolling-row DP using intermediate values. Our caching blends final profiles AFTER the DP. Mathematically equivalent but different float accumulation paths across 35 merge steps.
- **`match_calc` sparse representation**: C uses `cpmxpd`/`cpmxpdn` sparse arrays to skip zero-frequency residues during score computation. Our `match_score` iterates all alphabet positions. Same result but different float ordering.


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
