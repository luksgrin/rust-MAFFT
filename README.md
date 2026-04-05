# MAFFT-rs

A Rust reimplementation of [MAFFT](https://mafft.cbrc.jp/alignment/software/), the widely-used multiple sequence alignment tool originally written in C by Kazutaka Katoh.

This project provides a single `mafft-rs` binary that reads FASTA sequences and produces a multiple sequence alignment, replacing the original C toolchain (5+ binaries orchestrated by a 3000-line shell script) with a standalone Rust executable.

## Status

**Working prototype.** The core alignment pipeline (progressive alignment, iterative refinement, FFT-accelerated homology detection) is implemented and produces valid alignments for protein and DNA sequences. Alignment quality is approximately 50% of the original C implementation's sum-of-pairs score on the included test dataset — see [Known Limitations](#known-limitations) for details.

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
# Fast unit tests (all crates)
cargo test --workspace --lib

# Integration tests (requires release build, ~5s)
cargo test -p mafft-core --release --test end_to_end

# Build and test C reference (optional)
make -C mafft-upstream/core
```

## Usage

```bash
# Basic usage (FFT-NS-2, fast default)
mafft-rs sequences.fasta > aligned.fasta

# Read from stdin
cat sequences.fasta | mafft-rs > aligned.fasta

# Write to file
mafft-rs -o aligned.fasta sequences.fasta

# Quiet mode (suppress progress messages)
mafft-rs -q sequences.fasta > aligned.fasta
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
| `mafft --clustalout input.fa` | `mafft-rs --format clustal input.fa` | Clustal output |
| `mafft --phylipout input.fa` | `mafft-rs --format phylip input.fa` | PHYLIP output |

### Features not yet supported

| Original MAFFT flag | Status |
|---------------------|--------|
| `--add`, `--addfragments` | Not implemented (adding sequences to existing alignment) |
| `--parttree`, `--dpparttree` | Not implemented (PartTree for 10K+ sequences) |
| `--allowshift` | Not implemented (warp/shift gap penalty) |
| `--thread N` | Not implemented (multithreading) |
| `--nofft` | Not needed (FFT is used automatically when beneficial) |
| `--retree N` | Accepted but not yet effective |
| `--op`, `--ep`, `--bl` | Penalty tuning not exposed |
| `--kimura N` | Distance model tuning not exposed |
| RNA modes (`--qinsi`, `--xinsi`) | Not implemented |
| Structure modes (`--scarnalike`) | Not implemented |

## Architecture

The project is organized as a Cargo workspace with 9 crates:

```
crates/
  mafft-sys/       Raw FFI bindings to MAFFT C code (for cross-validation)
  mafft-types/     Shared Rust types (HomologyRegion, Sequence, ScoringContext, etc.)
  mafft-io/        FASTA, Clustal, PHYLIP, hat2 I/O
  mafft-scoring/   Substitution matrices (BLOSUM, JTT, TM, DNA) and gap penalties
  mafft-fft/       FFT-based homology detection (using rustfft + num-complex)
  mafft-align/     Pairwise and profile alignment algorithms (NW, SW, generalized affine)
  mafft-tree/      Distance computation, NJ, UPGMA, guide tree construction
  mafft-core/      Progressive alignment engine, iterative refinement, MafftEngine
  mafft-bin/       CLI binary (mafft-rs)
```

### How it was built

The migration followed a bottom-up, phase-by-phase strategy:

1. **Phase 0**: Cargo workspace + FFI bindings to C code + shared Rust types
2. **Phase 1**: FASTA/Clustal/PHYLIP I/O (lenient parser for MAFFT's non-standard headers)
3. **Phase 2**: Scoring matrices (BLOSUM 30-80, JTT, TM, DNA with IUPAC ambiguity codes)
4. **Phase 3**: FFT engine (replaced hand-rolled Cooley-Tukey with `rustfft`, multi-channel per-residue-type correlation)
5. **Phase 4a**: Pairwise alignment (Needleman-Wunsch, Smith-Waterman, generalized affine gap)
6. **Phase 4b**: Tree construction (NJ, UPGMA, production MUSCLE-style builder with nearest-neighbor heuristic)
7. **Phase 5**: Profile alignment (group-to-group DP, FFT-accelerated anchoring, local homology constraints)
8. **Phase 6**: Progressive alignment engine + iterative refinement (TreeDependentIteration port)
9. **Phase 7**: CLI binary consolidating all C entry points into a single `mafft-rs`

Each phase was tested independently against the C reference implementation, maintaining a suite of 96 tests (92 unit + 4 integration).

### Key design decisions

- **No global state.** The C code uses ~400 `extern` globals. Rust modules use owned `ScoringContext`, `Topology`, `Profile` structs passed explicitly.
- **`num_complex::Complex64`** replaces the C `Fukusosuu` struct. `rustfft` replaces the hand-rolled Cooley-Tukey FFT.
- **Per-group gap stripping** strips globally-all-gap columns before profile alignment (conservative but correct; see TODO.md for the more aggressive per-group approach).
- **Three-way gap insertion** (matching C's `insertnewgaps()`): group1 follows cursor1, group2 follows cursor2, "other" sequences follow cursor1 with gaps at Insert positions.

## Known limitations

### Alignment quality

On the included 36-sequence protein test dataset (`mafft-upstream/test/sample`), the Rust implementation achieves approximately 52% of the C implementation's sum-of-pairs identity score. This gap is due to:

- **Conservative gap stripping**: we strip only globally-all-gap columns, while C strips per-group, producing shorter (faster, better) profiles.
- **Simplified FFT anchoring**: the anchor selection heuristic is less tuned than C's.
- **Scoring normalization differences**: subtle differences in how scoring matrices are scaled and applied.

The alignment is structurally correct (all sequences have the same width, ungapped sequences match originals, residue content is preserved).

### Performance

- Single-threaded only (the C version uses pthreads). `rayon` integration is planned.
- Per-group gap stripping not yet implemented (DP matrices can be larger than necessary).
- No SIMD optimization for inner DP loops.

### Missing features

See the [Features not yet supported](#features-not-yet-supported) table above and `TODO.md` for detailed technical notes.

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

## License

- Rust code (`crates/`): [MIT](LICENSE-MIT)
- Original MAFFT C code (`mafft-upstream/`): [BSD-3-Clause](mafft-upstream/license)

## Credits

- Original MAFFT by [Kazutaka Katoh](https://mafft.cbrc.jp/alignment/software/) (CBRC, AIST)
- Rust reimplementation by Lucas Goiriz (luksgrin)
