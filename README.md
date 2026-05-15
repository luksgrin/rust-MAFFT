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
- L-INS-1, G-INS-1, E-INS-1 (`--localpair`/`--globalpair`/`--genafpair` with `--maxiterate 0`)
- L-INS-i, G-INS-i, E-INS-i (with iterative refinement)
- BLOSUM 30 / 45 / 50 / 62 / 80 (FFT and NW)
- JTT 100 / 200, TM 100 / 200 (FFT and NW)
- `--parttree`, `--dpparttree`, `--parttree --nofft` (divide-and-conquer guide tree)
- `--add`, `--add --nofft`, `--add --keeplength`
- Q-INS-i (RNA, requires `mxscarnamod`)
- `--allowshift --globalpair --maxiterate 0` (closed 2026-05-12)
- `--reorder` / `--inputorder` for all modes including PartTree (closed 2026-05-13)
- `--treeout` for all modes including `--parttree` and `--dpparttree` (closed 2026-05-13)

Every progressive merge step matches in score and width and every refinement iteration converges to C's exact alignment.

### Parity matrix (36-seq protein sample, 2026-05-12)

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
| TM 100 / 200 (NW + FFT)                     | match   | match      | 0    | ✓ byte-exact |
| `--parttree`, `--dpparttree`                | 752     | 752        | 0    | ✓ byte-exact |
| `--parttree --nofft`                        | 752     | 752        | 0    | ✓ byte-exact (closed 2026-05-12) |
| `--parttree --reorder`                      | 752     | 752        | 0    | ✓ byte-exact (closed 2026-05-13) |
| `--add` (30+6 fixture)                      | 741     | 741        | 0    | ✓ byte-exact |
| `--add --nofft` (30+6 fixture)              | 741     | 741        | 0    | ✓ byte-exact |
| `--add --keeplength` (30+6 fixture)         | 595     | 595        | 0    | ✓ byte-exact |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0†   | ✓ byte-exact (case-insensitive) |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 0†   | ✓ byte-exact (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0    | ✓ byte-exact (closed 2026-05-12) |
| `--reorder` (non-PartTree modes)            | match   | match      | 0    | ✓ byte-exact (closed 2026-05-13) |

† We uppercase residues; C preserves case. With `diff -i` (case-insensitive) RNA produces 0 lines.

Every mainstream mode is byte-identical to C MAFFT 7.526. The remaining
items in `TODO.md` are latent / coverage / performance gaps and missing
CLI flags, not active divergences.

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
# Full Rust suite (all crates, 269 tests as of 2026-05-12)
cargo test --workspace --exclude pymafft --release

# Just unit tests (faster)
cargo test --workspace --exclude pymafft --lib

# End-to-end integration tests (58 tests, requires release build)
cargo test -p mafft-core --release --test end_to_end

# Build C reference for cross-validation (optional)
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
| `--allowshift` | Supported (byte-identical to C, closed 2026-05-12) |
| `--nofft` | Supported |
| `--retree N` | Supported (default 2 for FFT-NS, 1 for INS-i, matching C) |
| `--op`, `--ep`, `--bl` | Supported |
| `--thread N` | Supported (Rayon) |
| `--kimura N` | Supported (custom Kimura R for DNA PAM generation) |
| `--parttree`, `--dpparttree` | Supported (divide-and-conquer tree for 10K+ sequences) |
| `--groupsize N` | Supported (partition size for PartTree, default 150) |
| `--qinsi` (Q-INS-i) | Supported (requires `mxscarnamod` in PATH) |
| `--xinsi` (X-INS-i) | Supported (requires `contrafold` in PATH) |
| `--scarnalike` | Supported (requires `dash_client` in PATH) |
| `--reorder` / `--inputorder` | Supported for non-PartTree modes (byte-identical to C MAFFT 7.526) |
| `--treeout` | Supported for all modes including `--parttree`/`--dpparttree` (byte-identical to C MAFFT 7.526) |
| `--auto`, `--seed`, `--treein`, `--memsave`, `--anysymbol`, `--leavegappyregion` | **Not yet implemented** (see `TODO.md` §B.3) |

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

Current counts as of 2026-05-15 (`cargo test --workspace --exclude pymafft --release`: **323 passed, 0 failed, 0 ignored**):

| Suite | Count | What |
|-------|-------|------|
| Rust unit tests (`--lib`) | 150 | All crates, all modules (incl. 3 new `msalign` Hirschberg DP tests) |
| Rust binary tests (`mafft-rs`) | 10 | CLI helper functions (`decide_auto`, `replace_unusual`, etc.) |
| Rust integration tests (`end_to_end`) | 77 | Byte-level parity with C across all supported modes + DP diagnostics |
| Rust FFI cross-validation tests | 61 | Cell-by-cell matrix equality, single-pair `G__align11` / `A__align` / `genL__align11` / warp DP / FFT equality, PartTree pipeline (8 tests), memsavetree (4 tests), MSalignmm profile alignment |
| Rust other integration tests | 25 | `trace_refinement` (15), `integration` mafft-io (7), `imp_*` (3) |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total Rust** | **323** | |

Regression guards for C parity are in `crates/mafft-core/tests/end_to_end.rs` (mode-level byte-identity), `crates/mafft-core/tests/cross_validate_*.rs` (FFI-level cell/function equality), and `crates/mafft-tree/tests/cross_validate_parttree.rs` (PartTree pipeline equality). Any regression in DP indexing, boundary handling, FFT anchor segment gaps, retree distance, refinement-tree distance, or pairwise/profile consistency will fail at least one of these.

## Known limitations

### `--xinsi` (untestable)

Wired correctly but requires Stanford's `CONTRAfold v2.02+` binary, which is not shipped by upstream MAFFT. The binary emits a "contrafold not found" diagnostic when run without the external dependency. See `TODO.md` §D.

### Missing CLI flags

Only `--seedtable` remains unimplemented — it produces an "unknown argument" error. See `TODO.md` §B.3.

`--memsave` and `--nomemsave` are accepted CLI shims (with C MAFFT's gating against `--localpair`/`--globalpair`/`--genafpair`/etc.). For inputs that fit in memory (≤ 30000 in length per sequence) the alignment is byte-identical to C MAFFT's `--memsave` output. The Hirschberg-style linear-space DP (`mafft_align::msalignmm`) is ported and verified to return optimal scores on identical-input and gap-required-with-unique-optimum tests, but isn't yet wired into the engine because matching C's exact tie-break behavior on score-tied alignments needs FFI cross-validation. See `TODO.md` §B.3 for the residual.

Closed since 2026-05-13: `--reorder`/`--inputorder` (§AA), `--treeout` (§AB), `--treein` (§AC), `--auto` (§AD), `--memsavetree` (§AE), `--anysymbol`/`--preservecase` (§AF), `--leavegappyregion`/`--legacygappenalty` (§AG), `--seed` (§AH), and `--memsave`/`--nomemsave` (§B.3) — all byte-identical to C MAFFT 7.526. PartTree's two-pass reorder (`splittbfast` CALL 1 with raw 6-mer distances, CALL 2 with `naivepairscore11` on the first-pass alignment) is composed as `final_order[k] = call1_order[call2_order[k]]`.

### Case preservation

C preserves the input case of residues (e.g., lowercase RNA). We uppercase all residues before alignment. The alignment itself (gap placement) is identical — `rna_nofft_case_insensitive_identical_to_c` verifies this.

### Performance

- Per-group gap stripping in progressive alignment is not implemented; we strip globally instead. Functionally correct but wastes DP work. See `TODO.md` §C.1.
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
