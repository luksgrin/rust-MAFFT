# MAFFT-rs

A Rust reimplementation of [MAFFT](https://mafft.cbrc.jp/alignment/software/), the widely-used multiple sequence alignment tool originally written in C by Kazutaka Katoh.

This project provides:
- **`mafft-rs`** — a single CLI binary replacing the original 5+ C binaries and shell script
- **`mafft-core`** — a Rust library crate for programmatic use
- **`pymafft`** — Python bindings via PyO3

## Status

**Working implementation.** The core alignment pipeline (progressive alignment, iterative refinement, FFT-accelerated homology detection) is implemented and produces valid alignments for protein and DNA sequences. On the included 36-sequence protein test dataset, the following mode/flag combinations produce **byte-identical output to C MAFFT 7.526** (`diff rust.fa c.fa` returns 0 lines):

- FFT-NS-2 (default), NW-NS-2 (`--nofft`)
- FFT-NS-i (`--maxiterate 100`)
- L-INS-1 (`--localpair --maxiterate 0`) — closed 2026-05-05
- BLOSUM 30 / 45 / 62 / 80 (FFT and NW)
- JTT 200 (FFT)
- TM 100 / 200 (NW)
- RNA NW (case-insensitive, since we uppercase input)

Every progressive merge step matches in score and width and every refinement iteration converges to C's exact alignment for the byte-exact modes. See [Known limitations](#known-limitations) and `TODO.md` for the remaining mode/flag combinations that still diverge.

### Parity matrix (36-seq protein sample)

| Mode / flags                                | C width | Rust width | diff | Status |
|---------------------------------------------|---------|------------|------|--------|
| FFT-NS-2 (default)                          | 717     | 717        | 0    | ✓ byte-exact |
| NW-NS-2 (`--nofft`)                         | 717     | 717        | 0    | ✓ byte-exact |
| FFT-NS-i (`--maxiterate 100`)               | 721     | 721        | 0    | ✓ byte-exact |
| L-INS-1 (`--localpair --maxiterate 0`)      | 719     | 719        | 0    | ✓ byte-exact (closed 2026-05-05) |
| L-INS-i (`--localpair`)                     | 719     | 725        | 934  | ✗ refinement loop diverges |
| G-INS-1 (`--globalpair --maxiterate 0`)     | 746     | 747        | 974  | ✗ progressive diverges |
| G-INS-i (`--globalpair`)                    | 746     | 714        | 939  | ✗ progressive + refinement |
| E-INS-1 (`--genafpair --maxiterate 0`)      | 729     | 719        | 948  | ✗ uses local pairwise instead of generalized-affine |
| E-INS-i (`--genafpair`)                     | 729     | 725        | 984  | ✗ same as above + refinement |
| BLOSUM 30 / 45 / 62 / 80 (FFT)              | match   | match      | 0    | ✓ byte-exact |
| BLOSUM 80 NW (`--bl 80 --nofft`)            | 712     | 712        | 0    | ✓ byte-exact |
| BLOSUM 50 FFT                               | 712     | 738        | 961  | ✗ FFT multi-lag tie-break |
| JTT 200 (FFT)                               | 729     | 729        | 0    | ✓ byte-exact |
| JTT 100 (FFT)                               | 732     | 732        | 4    | ✗ FFT tie-break (same column shift in 1 seq) |
| TM 200 NW                                   | 765     | 765        | 0    | ✓ byte-exact |
| TM 100 NW                                   | 767     | 767        | 0    | ✓ byte-exact |
| TM 200 FFT                                  | 767     | 765        | 148  | ✗ FFT tie-break |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0†   | ✓ byte-exact (case-insensitive) |
| `--parttree --nofft`                        | —       | —          | ~940 | ✗ subtree-grouping bug |
| `--add --nofft`                             | —       | —          | ~900 | ✗ insertnewgaps incomplete |

† We uppercase residues; C preserves case. With `diff -i` the RNA path produces 0 lines.

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
# Full Rust suite (all crates, 223 tests)
cargo test --workspace --exclude pymafft --release

# Just unit tests (faster)
cargo test --workspace --exclude pymafft --lib

# End-to-end integration tests (39 tests, requires release build)
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

Current counts as of 2026-05-05 (`cargo test --workspace --exclude pymafft --release`: 223 passed, 0 failed, 0 ignored):

| Suite | Count | What |
|-------|-------|------|
| Rust unit tests | ~140 | All crates, all modules (per-crate `--lib` runs) |
| Rust integration tests (`end_to_end`) | 39 | Byte-level parity with C (FFT-NS-2, NW-NS-2, FFT-NS-i, `--bl 30/45/80` with/without FFT, `--jtt 200`, `--tm 100/200 --nofft`, RNA `--nofft`), DP diagnostics |
| Rust FFI cross-validation tests | 16 | Cell-by-cell matrix equality and constrained-DP equivalence vs C via FFI (BLOSUM45/50/62/80, JTT 200, JTT 100, TM 200 n_dis + n_disFFT, BL50 n_disFFT, DNA, plus 5 constrained-align tests in `cross_validate_constrained_align.rs`) |
| C alignment tests | 8 | FFT-NS-2, FFT-NS-i, G-INS-i, L-INS-i, parttree, etc. |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total Rust** | **223** | |

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

### Iterative refinement (FFT-NS-i, L-INS-i, G-INS-i, E-INS-i)

FFT-NS-i (`--maxiterate 100`) is byte-identical to C — see `fftnsi_byte_identical_to_c`. The fix that closed the last 1-column gap was that the mafft script's second `dndpre` invocation (which writes the hat2 file `dvtditr` reads) is called WITHOUT `-h 0`, so it uses `poffset = -123` → matrix shifted by `+73` versus the matrix DP uses. Our refinement branch in `engine.rs` builds a `+73`-shifted copy of the substitution matrix for the refinement-tree distance computation only.

**L-INS-i without refinement (`--localpair --maxiterate 0`) is byte-identical to C as of 2026-05-05.** The closure required three fixes:
1. `opt` field of homology regions was stored on the pre-rescale (hat3-file) side of C's `tbfast.c:2202` rescale by `* 600 / 5.8`, producing impmtx values 100× smaller than C's. Now stored as `isumscore / sumoverlap` (post-rescale value) directly.
2. The L-INS-i path was rounding distances to 3 decimals (mimicking the hat2 `%#6.3f` round-trip). That's correct for FFT-NS-i + dndpre but wrong for L-INS-i: `tbfast` with `callpairlocalalign=1` keeps full-precision `iscore[]` in memory. Removed the rounding for the constraint-aware path.
3. The CLI was treating `--maxiterate 0` as "use the mode default 1000" because `args.maxiterate` was a `usize` defaulting to 0 and could not distinguish "unset" from "explicit 0". Switched to `Option<usize>`.

**L-INS-i with refinement (`--localpair`, default 1000 iterations), G-INS-i, E-INS-i still diverge** — the refinement loop with constraints picks different rearrangements than C's. See `TODO.md` §3 for the diagnostic plan. G-INS-1 (`--globalpair --maxiterate 0`) and E-INS-1 (`--genafpair --maxiterate 0`) also diverge even without refinement: G-INS-1 likely needs the same opt/distance fixes applied to the global-pair path; E-INS-i still uses Smith-Waterman pairwise (it should use generalized-affine — see TODO §2).

### Non-default BLOSUM matrices (`--bl N`)

`--bl 30`, `--bl 45`, `--bl 62` (default), `--bl 80` are all byte-identical to C — see `fftns2_bl30_byte_identical_to_c`, `fftns2_bl45_byte_identical_to_c`, `fftns2_bl80_byte_identical_to_c`, `nofft_bl80_byte_identical_to_c`. MAFFT's BLOSUM45 and BLOSUM80 tables differ from standard NCBI at a handful of cells (BL45: P/H = -2; BL80: H/R = 0, F/M = 0, P/R = -3, V/I = 4) — matching upstream MAFFT byte-for-byte requires keeping those variant cells. `--bl 50`'s 26×26 substitution matrix matches C cell-by-cell (`cross_validate_blosum50_n_dis_cell_by_cell`, `cross_validate_blosum50_n_dis_fft_cell_by_cell`) but its end-to-end alignment diverges from C (Rust width 738 vs C 712 on the 36-seq test); the DP is selecting a different equivalent path because our `find_fft_anchors` picks the single best lag while C's `Falign` accumulates segments from all `NKOUHO=20` candidate lags. See `TODO.md` §4.

### Substitution-model flags (`--jtt`, `--tm`)

`--jtt N` and `--tm N` are wired through the CLI (matching `mafft --jtt N` / `mafft --tm N`). The PAM number selects how many iterations of the JTT one-step transition matrix are multiplied together before the log-odds transform.

- `--jtt 200` is byte-identical to C in FFT-NS-2 (`fftns2_jtt200_byte_identical_to_c`).
- `--tm 100 --nofft` and `--tm 200 --nofft` are byte-identical to C (`nofft_tm100_byte_identical_to_c`, `nofft_tm200_byte_identical_to_c`).
- `--tm` with FFT and `--jtt 100 (FFT)` produce alignments that score identically but place gaps differently in a few positions. Cell-by-cell matrix matches C exactly (FFI cross-validation tests `cross_validate_jtt_n_dis`, `cross_validate_jtt100_n_dis`, `cross_validate_tm_n_dis`, `cross_validate_tm_n_dis_fft`); the divergence is the same FFT multi-lag accumulation issue as `--bl 50`. See `TODO.md` §4-§5.

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
