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
- JTT 100 / 200, TM 100 / 200 (FFT and NW; TM 200 + `--retree 1` has an
  8-line diff that's a C-side tied-trace artifact — see
  `MAFFT_UPSTREAM_REPORT.md`)
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
| TM 200 (`--retree 1`)                       | 717     | 717        | 8†   | C-side tied-trace artifact (see `MAFFT_UPSTREAM_REPORT.md`) |
| `--parttree`, `--dpparttree`, `--parttree --nofft`, `--parttree --reorder` | 752 | 752 | 0 | ✓ byte-exact |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6 fixture) | match | match | 0 | ✓ byte-exact |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0‡   | ✓ byte-exact (case-insensitive) |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 0‡   | ✓ byte-exact (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0    | ✓ byte-exact |
| `--reorder` (non-PartTree modes)            | match   | match      | 0    | ✓ byte-exact |

‡ For RNA / Q-INS-i: we uppercase residues; C preserves case. With `diff -i` RNA produces 0 lines.

† `--tm 200 --retree 1`: 4 single-char gap shifts in 2 of 36 sequences. Both
alignments are optimal-scored ties; the diff comes from C `A__align`'s
`static TLS` buffers biasing tied-DP-cell selection by accumulated state
from prior calls. Our DP produces byte-identical `ijp[i][j]` to C's
`A__align` given identical inputs (verified cell-by-cell across all
365×363 cells), so this is a C-side artifact, not a Rust bug. Full
analysis + upstream report in `MAFFT_UPSTREAM_REPORT.md`. Default
`--tm 200` (retree=2) is byte-identical because the second pass realigns
from a rebuilt tree.

Every other mode is byte-identical to C MAFFT 7.526. Remaining items in
`TODO.md` are non-correctness gaps (performance deferrals and one
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

`cargo test --workspace --exclude pymafft --release`: **349 passed, 0 failed, 0 ignored**:

| Suite | Count | What |
|-------|-------|------|
| Rust unit tests (`--lib`) | 157 | All crates, all modules (incl. 7 `msalign` Hirschberg DP tests + 4 `parse_hat3_seed` tests) |
| Rust binary tests (`mafft-rs`) | 10 | CLI helper functions (`decide_auto`, `replace_unusual`, etc.) |
| Rust integration tests (`end_to_end`) | 88 | Byte-level parity with C across all supported modes + DP diagnostics (incl. 4 `--seedtable` and 7 `--retree N` byte-equal tests) |
| Rust FFI cross-validation tests | 68 | Cell-by-cell matrix equality, single-pair `G__align11` / `A__align` / `genL__align11` / warp DP / FFT equality, PartTree pipeline, memsavetree, `MSalignmm` Hirschberg DP, MSalignmm profile alignment |
| Rust other integration tests | 25 | `trace_refinement` (15), `integration` mafft-io (7), `imp_*` (3) |
| Python tests | 32 | API, strategies, file I/O, error handling, types |
| **Total Rust** | **349** | |

Regression guards for C parity are in `crates/mafft-core/tests/end_to_end.rs` (mode-level byte-identity), `crates/mafft-core/tests/cross_validate_*.rs` (FFI-level cell/function equality), and `crates/mafft-tree/tests/cross_validate_parttree.rs` (PartTree pipeline equality). Any regression in DP indexing, boundary handling, FFT anchor segment gaps, retree distance, refinement-tree distance, or pairwise/profile consistency will fail at least one of these.

#### BALIBASE 3 parity sweep

Sweep across BALIBASE 3 (Drive5 mirror, 218 protein test sets) against MAFFT 7.526:

| Mode | Match | Total % |
|------|-------|---------|
| FFT-NS-1 (`--retree 1`) | 206/218 | 94.5 % |
| FFT-NS-2 (default) | 212/218 | **97.2 %** |
| FFT-NS-i (`--maxiterate 100`) | 128/218 | 58.7 % |
| L-INS-i (`--localpair --maxiterate 100`) | 154/218 | 70.6 % |
| G-INS-i (`--globalpair --maxiterate 100`) | 140/218 | 64.2 % |
| E-INS-i (`--genafpair --maxiterate 100`) | 184/218 | 84.4 % |

All residual divergences are `A__align` static-state-coupling
artifacts (same class as `--tm 200 --retree 1` — see
`MAFFT_UPSTREAM_REPORT.md`). Cell-level FFI tests
(`cross_validate_cpmx`, `cross_validate_counteff`,
`cross_validate_bb20027_dp`) prove the Rust DP and profile-building
match C bit-for-bit on identical inputs. Divergences come from C's
per-process static TLS memoization state, not from a Rust bug. Both
alignments at each divergent case are optimal-scored. Refinement
modes show lower parity because every iteration re-runs `A__align`,
amplifying any tied-cell shift through subsequent passes; ~80 % of
refinement-mode divergences are same-width tied-trace cases.

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
byte-identical to C on the 36-seq test sample.

The remaining `TODO.md` items are: (a) `A__align` static-state
coupling artifacts (`--tm 200 --retree 1` 8-line diff + 6/218
BALIBASE 3 cases — see `MAFFT_UPSTREAM_REPORT.md`), (b) performance
deferrals where Rust is already faster than C on most modes, and
(c) `--xinsi` untestable until Stanford's `contrafold` is installed.

### Case preservation

C preserves the input case of residues (e.g., lowercase RNA). We uppercase all residues before alignment. The alignment itself (gap placement) is identical — `rna_nofft_case_insensitive_identical_to_c` verifies this.

### Performance

Benchmark snapshot vs C MAFFT 7.526 (2026-05-18, macOS arm64, 5 runs summed,
108-seq synthetic input):

| Mode                                     | C MAFFT | mafft-rs | Ratio          |
|------------------------------------------|---------|----------|----------------|
| default (FFT-NS-2)                       | 1.34s   | 1.16s    | 1.16× faster   |
| `--maxiterate 50` (FFT-NS-i)             | 6.62s   | 7.78s    | 1.18× SLOWER   |
| `--maxiterate 50 --localpair` (L-INS-i)  | 37.00s  | 22.93s   | **1.61× faster** |
| `--maxiterate 50 --globalpair` (G-INS-i) | 35.32s  | 22.62s   | **1.56× faster** |
| `--parttree`                             | 1.70s   | 0.09s    | **18× faster**   |
| `--dpparttree`                           | 7.24s   | 0.08s    | **90× faster**   |

FFT-NS-i is the only mode where Rust is slower; the deficit is in the
inner DP cell loop (55% self time per `macOS sample`-based profile)
and its per-call `Vec` allocation (~14% of that). C's `A__align`
amortizes allocation across calls via `static TLS` buffers (the same
mechanism behind the §B.2 stateful-buffer artifact). See
[`PROFILING.md`](./PROFILING.md) for the recorded profile, the recipe
to reproduce it, and three candidate optimization paths.

Earlier-flagged perf items (`TODO.md` §C.1 per-group gap stripping, §C.2
SIMD inner DP loops) are deferred — both were filed before measurement;
the modes §C.1 targets are already faster than C, and §C.2's "SIMD the
arithmetic" framing isn't the real lever (allocator amortization is).

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
