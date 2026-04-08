# MAFFT-rs

A Rust reimplementation of [MAFFT](https://mafft.cbrc.jp/alignment/software/), the widely-used multiple sequence alignment tool originally written in C by Kazutaka Katoh.

This project provides:
- **`mafft-rs`** — a single CLI binary replacing the original 5+ C binaries and shell script
- **`mafft-core`** — a Rust library crate for programmatic use
- **`pymafft`** — Python bindings via PyO3

## Status

**Working prototype.** The core alignment pipeline (progressive alignment, iterative refinement, FFT-accelerated homology detection) is implemented and produces valid alignments for protein and DNA sequences. Alignment quality is approximately 62% of the original C implementation's sum-of-pairs score on the included test dataset — see [Known Limitations](#known-limitations) and [Gap Analysis](#gap-analysis-vs-original-mafft) for details.

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
| `mafft --nofft input.fa` | `mafft-rs --nofft input.fa` | Disable FFT (NW-NS-2) |
| `mafft --clustalout input.fa` | `mafft-rs --format clustal input.fa` | Clustal output |
| `mafft --phylipout input.fa` | `mafft-rs --format phylip input.fa` | PHYLIP output |

### Features not yet supported

| Original MAFFT flag | Status |
|---------------------|--------|
| `--add`, `--addfragments` | Not implemented (adding sequences to existing alignment) |
| `--parttree`, `--dpparttree` | Not implemented (PartTree for 10K+ sequences) |
| `--allowshift` | Not implemented (warp/shift gap penalty) |
| `--nofft` | Supported (disables FFT, forces pure DP) |
| `--retree N` | Supported (default 2, matching C) |
| `--op`, `--ep`, `--bl` | Supported |
| `--kimura N` | Accepted (stored, pending DNA PAM generation) |
| RNA modes (`--qinsi`, `--xinsi`) | Not implemented |
| Structure modes (`--scarnalike`) | Not implemented |

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
| Rust unit tests | 91 | All crates, all modules |
| Rust integration tests | 4 | End-to-end on real data + C reference comparison |
| C alignment tests | 8 | FFT-NS-2, FFT-NS-i, G-INS-i, L-INS-i, parttree, etc. |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total** | **136** | |

## Known limitations

### Alignment quality

On the included 36-sequence protein test dataset (`mafft-upstream/test/sample`), the Rust implementation achieves approximately 62% of the C implementation's sum-of-pairs identity score. The alignment is structurally correct (all sequences have the same width, ungapped sequences match originals, residue content is preserved). See the gap analysis below for the specific causes and fixes.

### Performance

- Per-group gap stripping not yet implemented (DP matrices can be larger than necessary). See `TODO.md`.
- No SIMD for the DP fill loops themselves (data dependencies prevent vectorization without anti-diagonal restructuring).

## Gap analysis vs original MAFFT

The following items explain the quality and feature gap between `mafft-rs` and the original C MAFFT, ordered by impact on alignment quality:

### High impact (alignment quality)

| Item | Description | Effort |
|------|-------------|--------|
| **Per-group gap stripping** | C strips columns all-gap within each group independently before profile alignment, producing shorter profiles. We strip only globally-all-gap columns. This directly affects profile quality and DP cost. | Medium — algorithm is correct but re-insertion interleaving has a bug on large inputs (see `TODO.md`) |
| **FFT anchor quality** | The C code's `seq_vec_3` vectorization, `getKouho` candidate selection, and `blockAlign2` anchor pairing are more tuned. Our multi-channel FFT works but anchor placement may differ. | Medium |
| **Guide tree fidelity** | Our `musclesupg` is equivalent algorithmically but may produce different trees due to tie-breaking or floating-point order differences. | Low — compare tree topologies |
| **`commongappick` during refinement** | C strips common gaps before each refinement re-alignment. Our refinement uses the full-width profiles. | Medium — same per-group stripping issue |
| **Distance computation for guide tree** | C uses 6-tuple distance with memoized frequency tables. Our `ktuple_distance` uses HashMap, which may produce slightly different values. | Low |

### Medium impact (missing features)

| Item | Description | Effort |
|------|-------------|--------|
| **`--add` / `--addfragments`** | Adding new sequences to an existing alignment (`addonetip` in C). Used frequently in incremental workflows. | Medium |
| **`--retree N`** | Rebuilding the guide tree N times. Currently accepted but ignored (always builds once). | Low |
| **`--op`, `--ep`, `--bl`, `--kimura`** | Gap penalty and distance model tuning. The parameters exist internally but aren't exposed via CLI. | Low |
| **`n_disLN` matrix** | Log-normal scoring variant used in a specific code path. | Low |
| **`ribosumdis[37][37]`** | RNA ribosum composite matrix. Raw data present but 37x37 assembly not done. | Low |
| **Warp/shift (`--allowshift`)** | Extra DP state for long-range gap shifts. | Medium |

### Low impact (niche features)

| Item | Description | Effort |
|------|-------------|--------|
| **PartTree (`--parttree`, `--dpparttree`)** | Divide-and-conquer tree for 10K+ sequences. | High |
| **RNA modes (`--qinsi`, `--xinsi`)** | McCaskill/CONTRAfold RNA structure integration. | High — requires external tools |
| **Structure alignment (`--scarnalike`)** | 3D structure-aware alignment via DASH. | High — requires external tools |
| **`veryfastsupg_int`** | Fast integer-distance UPGMA variant. | Low |
| **`blockAlign3`** | O(n^2) anchor selection variant. | Low |

### Path to 100% replication

To achieve byte-for-byte identical output with the C implementation on all test cases, three categories of work are needed: quality fixes (close the SP score gap), behavioral parity (match C's exact logic), and feature completeness (support all C flags).

#### Quality fixes (close the 52% → 100% SP score gap)

| # | Item | Description | Effort |
|---|------|-------------|--------|
| 1 | **Fix per-group gap stripping** | The single highest-impact item. C strips columns all-gap within each group independently (`commongappick`), producing shorter profiles and better DP. Our global stripping produces profiles 2-5x longer. The algorithm is correct but the re-insertion interleaving has a bug on large inputs (see `TODO.md`). | Medium |
| 2 | **Guide tree topology comparison** | Run both C and Rust on `test/sample`, compare join order and branch lengths. A different tree produces a fundamentally different progressive alignment. | Low |
| 3 | **FFT anchor position comparison** | For each progressive merge step, compare the anchor positions chosen by C vs Rust. Wrong anchors = wrong segment boundaries = wrong sub-alignments. | Medium |

#### Behavioral parity (match C's exact output)

| # | Item | Description | Effort |
|---|------|-------------|--------|
| 4 | **`commongappick` during refinement** | C strips common gaps before each refinement re-alignment. Our refinement uses full-width profiles. Same per-group stripping fix as item 1. | Medium |
| 5 | **`n_disLN` matrix** | Log-normal scoring variant used in a specific refinement code path. | Low |
| 6 | **Distance matrix: match C's 6-tuple counting** | C uses memoized frequency tables with `commonsextet_p`. Our HashMap-based k-tuple may produce slightly different values. | Low |
| 7 | **Match C's output ordering and formatting** | C outputs sequences in input order with specific line wrapping and name formatting. Small formatting differences exist. | Low |

#### Feature completeness (support all C modes and flags)

| # | Item | Description | Effort |
|---|------|-------------|--------|
| 8 | **`--add` / `--addfragments`** | Add new sequences to an existing alignment (`addonetip` algorithm). Frequently used in incremental workflows. | Medium |
| 9 | **`--allowshift`** | Warp/shift gap penalty — extra DP state for long-range jumps (`penalty_shift`). | Medium |
| 10 | **`--parttree` / `--dpparttree`** | PartTree divide-and-conquer for 10K+ sequence datasets. | High |
| 11 | **`ribosumdis[37][37]`** | Assemble the 37×37 RNA ribosum composite matrix from the 4×4 and 16×16 components already in `dna.rs`. | Low |
| 13 | **RNA modes (`--qinsi`, `--xinsi`)** | Integrate McCaskill/CONTRAfold RNA secondary structure predictions into alignment scoring. | High |
| 14 | **Structure alignment (`--scarnalike`)** | 3D structure-aware alignment via DASH client. | High |
| 15 | **`veryfastsupg_int`** | Fast integer-distance UPGMA variant (performance optimization). | Low |
| 16 | **`blockAlign3`** | O(n²) anchor selection variant (rarely triggered). | Low |

#### Summary

- **Items 1-3**: Close the quality gap (currently 47% SP ratio — temporary dip as individual pieces are fixed to match C exactly; quality will converge as more pieces are corrected).
- **Items 4-7**: Achieve behavioral parity (exact output matching on standard benchmarks).
- **Items 8-15**: Full feature completeness (all C flags supported).

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
