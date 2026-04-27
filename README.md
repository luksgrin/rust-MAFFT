# MAFFT-rs

A Rust reimplementation of [MAFFT](https://mafft.cbrc.jp/alignment/software/), the widely-used multiple sequence alignment tool originally written in C by Kazutaka Katoh.

This project provides:
- **`mafft-rs`** — a single CLI binary replacing the original 5+ C binaries and shell script
- **`mafft-core`** — a Rust library crate for programmatic use
- **`pymafft`** — Python bindings via PyO3

## Status

**Working implementation.** The core alignment pipeline (progressive alignment, iterative refinement, FFT-accelerated homology detection) is implemented and produces valid alignments for protein and DNA sequences. On the included 36-sequence protein test dataset, **FFT-NS-2 (default), NW-NS-2 (`--nofft`), FFT-NS-i (`--maxiterate 100`), `--bl 80` (with and without `--nofft`), `--jtt 200`, and `--tm 100/200 --nofft` all produce byte-identical output to C MAFFT 7.526** — every progressive merge step matches in score and width, every refinement iteration converges to C's exact alignment, and `diff rust_output.fa c_output.fa` returns 0 lines for all of these mode/scoring combinations. See [Known limitations](#known-limitations) for remaining gaps.

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
# Fast unit tests (all crates, 111 tests)
cargo test --workspace --exclude pymafft --lib

# Integration tests (requires release build, 26 tests)
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
| Rust unit tests | 111 | All crates, all modules |
| Rust integration tests | 37 | End-to-end on real data, byte-level parity with C (FFT-NS-2, NW-NS-2, FFT-NS-i, `--bl 80` with/without FFT, `--jtt 200`, `--tm 100/200 --nofft`), DP diagnostics |
| Rust FFI cross-validation tests | 8 | Cell-by-cell matrix equality vs C via FFI (BLOSUM62, BLOSUM80, JTT 200, JTT 100, TM 200 n_dis + n_disFFT, DNA) |
| C alignment tests | 8 | FFT-NS-2, FFT-NS-i, G-INS-i, L-INS-i, parttree, etc. |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total** | **190** | |

Regression guards for C parity (all in `crates/mafft-core/tests/end_to_end.rs`):

- `nofft_byte_identical_to_c` — final alignment matches C MAFFT 7.526's `--nofft` output byte-for-byte (fixture: `tests/fixtures/sample.nwns2`).
- `nofft_per_step_matches_c` — every per-merge `(clus1, clus2, width, score)` tuple matches C's across both retree passes (fixture: `tests/fixtures/sample.nwns2.steps`, 70 merges).
- `fftns2_byte_identical_to_c` — FFT-NS-2 (default strategy) output matches C's `mafft-upstream/test/sample.fftns2` reference byte-for-byte, guarding the full pipeline: FFT anchoring, segment gap handling, inter-anchor DP, retree distance, and UPGMA.
- `fftnsi_byte_identical_to_c` — FFT-NS-i (`--maxiterate 100`) output matches C's `mafft-upstream/test/sample.fftnsi` reference byte-for-byte (asserts both width equality and per-sequence equality), guarding the iterative-refinement pipeline end-to-end including the dndpre offset-shift step the mafft script applies before dvtditr.
- `fftns2_bl80_byte_identical_to_c` — FFT-NS-2 with `--bl 80` matches C byte-for-byte (fixture: `tests/fixtures/sample.bl80.fftns2`), guarding MAFFT's variant of the BLOSUM80 substitution-matrix table.
- `nofft_bl80_byte_identical_to_c` — NW-NS-2 with `--bl 80 --nofft` matches C byte-for-byte (fixture: `tests/fixtures/sample.bl80.nwns2`), same matrix-table guard via the non-FFT path.
- `fftns2_jtt200_byte_identical_to_c` — FFT-NS-2 with `--jtt 200` matches C byte-for-byte (fixture: `tests/fixtures/sample.jtt200.fftns2`), guarding the JTT lower-triangle accepted-point-mutation table, the PAM matrix exponentiation loop, and the shared normalize/600-scale/offset pipeline.
- `nofft_tm200_byte_identical_to_c` — NW-NS-2 with `--tm 200 --nofft` matches C byte-for-byte (fixture: `tests/fixtures/sample.tm200.nwns2`), guarding the TM upper-triangle table in `tm_rsr_matrix` and the TM frequency vector. Before this guard, `--tm` silently produced JTT-like output because the upper triangle was never populated.
- `nofft_tm100_byte_identical_to_c` — `--tm 100 --nofft` exercises the same TM data at a non-default PAM, catching drift in the matrix-power path.

Plus three layers of unit-level guards in `crates/mafft-scoring/`:

- `blosum::blosum80_tests::blosum80_table_layout_and_diagonals` — pins all 20 diagonal entries of the raw 210-cell BLOSUM80 lower-triangle array.
- `blosum::blosum80_tests::blosum80_mafft_variant_cells` — pins the four cells (H/R, F/M, P/R, V/I) where MAFFT's `tmpmtx80` deliberately differs from standard NCBI BLOSUM80, so a future "fix" to the NCBI standard can't slip through silently.
- `cross_validate_blosum80_n_dis_cell_by_cell` — calls C's `constants()` via FFI with BLOSUM80 and asserts every cell of the 26×26 normalized `n_dis` matches our `substitution_matrix`, catching pipeline regressions in the average-subtract / 600-scale / offset-subtract chain that would otherwise affect every cell.
- `nofft_op_override_byte_identical_to_c` — NW-NS-2 with `--op 2.5` matches C byte-for-byte (fixture: `tests/fixtures/sample.nwns2.op25`), guarding the gap-opening override path.
- `nofft_ep_override_byte_identical_to_c` — NW-NS-2 with `--ep 0.5` matches C byte-for-byte (fixture: `tests/fixtures/sample.nwns2.ep05`), guarding the scoring-matrix offset override path.
- `rna_nofft_case_insensitive_identical_to_c` — RNA NW-NS-2 output matches C byte-for-byte after case normalization (fixture: `tests/fixtures/samplerna.nwns2`), guarding the nucleotide alignment path.
- `diagnostic_simple_offset` — the minimal `ACDE` vs `WWWWACDEWWWW` reproducer must produce the optimal `----ACDE----` alignment with 4 matches.
- `diagnostic_align11_vs_profile` — `pairwise_align11` and `profile_align` must agree on 1×1 inputs (same operations, same score).

Any regression in DP indexing, boundary handling, FFT anchor segment gaps, retree distance, refinement-tree distance, or pairwise/profile consistency will fail at least one of these.

## Known limitations

### Iterative refinement (FFT-NS-i, G-INS-i, L-INS-i, E-INS-i)

FFT-NS-i (`--maxiterate 100`) is now byte-identical to C — see `fftnsi_byte_identical_to_c`. The fix that closed the last 1-column gap was discovering that the mafft script's second `dndpre` invocation (which writes the hat2 file dvtditr reads) is called WITHOUT `-h 0`, so it uses `poffset = -123` → matrix shifted by `+73` versus the matrix DP uses. Our refinement branch in `engine.rs` builds a `+73`-shifted copy of the substitution matrix for the refinement-tree distance computation only.

`G-INS-i`, `L-INS-i`, and `E-INS-i` share the same refinement core but have not yet been cross-validated per-iteration against C. They likely benefit from the same fix; first canonical example to check is in TODO §5.

### Non-default BLOSUM matrices (`--bl N`)

`--bl 80` (with and without `--nofft`) is byte-identical to C — see `fftns2_bl80_byte_identical_to_c` and `nofft_bl80_byte_identical_to_c`. `--bl 30/45/50` share the same code path but don't yet have C-reference fixtures; trajectories are not formally validated.

### Substitution-model flags (`--jtt`, `--tm`)

`--jtt N` and `--tm N` are now wired through the CLI (matching `mafft --jtt N` / `mafft --tm N`). The PAM number selects how many iterations of the JTT one-step transition matrix are multiplied together before the log-odds transform.

- `--jtt 200` is byte-identical to C in FFT-NS-2 (`fftns2_jtt200_byte_identical_to_c`).
- `--tm 100 --nofft` and `--tm 200 --nofft` are byte-identical to C (`nofft_tm100_byte_identical_to_c`, `nofft_tm200_byte_identical_to_c`).
- `--tm` with FFT and `--jtt 100` produce alignments that score identically but place gaps differently in a few positions (DP tie-breaking). Cell-by-cell matrix matches C exactly in all four cases (FFI cross-validation tests `cross_validate_jtt_n_dis`, `cross_validate_jtt100_n_dis`, `cross_validate_tm_n_dis`, `cross_validate_tm_n_dis_fft`).

### PartTree (`--parttree`, `--dpparttree`)

`--parttree --nofft` diverges from C (~940-line diff). Needs investigation in `crates/mafft-tree/src/parttree.rs`.

### Adding sequences (`--add`, `--addfragments`, `--keeplength`)

`--add --nofft` diverges from C (~900-line diff). The add pipeline in `crates/mafft-core/src/add.rs` needs investigation against C's `addsingle`.

### Per-group gap stripping (performance)

The progressive alignment strips only columns that are all-gap across **all sequences globally**, instead of stripping per-group as C does. This is functionally correct (produces the same alignment) but wastes DP computation. See `TODO.md` for a detailed diagnosis and attempted approaches. The correct fix requires porting C's `insertnewgaps()` from `addfunctions.c`.

### Case preservation

C preserves the input case of residues (e.g., lowercase RNA). We uppercase all residues before alignment. The alignment itself (gap placement) is identical. The `rna_nofft_case_insensitive_identical_to_c` test verifies this.

### Performance

- No SIMD for the DP fill loops themselves (data dependencies prevent vectorization without anti-diagonal restructuring).


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
