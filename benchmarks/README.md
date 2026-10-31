# Benchmarks

Drivers that compare `mafft-rs` with C MAFFT on peak memory, thread
scaling and determinism, and that exercise `pymafft`'s `threads=` keyword.
They write tab-separated results to `benchmarks/results/`; the summaries
quoted in the documentation (README "Performance", `PROFILING.md`,
`docs/python-api/index.md`) are taken from the TSVs committed there.

| Script | Measures | Output |
|--------|----------|--------|
| `memory.py` | peak resident set size, C vs `mafft-rs`, per mode and `--thread N`, plus `--memsave` | `results/memory_<date>.tsv` |
| `threads.py` | determinism of `mafft-rs` across `--thread N`, C's run-to-run variation, wall-clock speedup | `results/threads_<date>.tsv` |
| `pymafft_threads.py` | `pymafft.align(..., threads=N)`: equality with the CLI, timing, GIL release | `results/pymafft_threads_<date>.tsv` |

`common.py` holds the shared paths, input list and the `/usr/bin/time`
parser; `summarize_memory.py` turns `memory.py` TSVs into the per-cell
Markdown table (`results/memory_<date>_summary.md`). For alignment-accuracy scores against the BAliBASE references see
`balibase_scores.py`.

## Reference build

Both sides are built from this checkout; no downloaded binary is used.

```sh
make -C mafft-upstream/core            # C reference -> mafft-upstream/binaries/
cargo build --release -p mafft-rs      # -> target/release/mafft-rs
```

C is run through the upstream wrapper `mafft-upstream/scripts/mafft` with
`MAFFT_BINARIES=mafft-upstream/binaries`, as in
`.github/workflows/build-test.yml`; the drivers set that variable
themselves. Both builds use the host's default floating-point contraction
policy (fused multiply-add on arm64, none on x86-64; see
`docs/architecture/byte-identity.md`), recorded in the `fp_policy` column.

For the Python driver:

```sh
python3 -m venv .venv
.venv/bin/pip install maturin
cd crates/pymafft && VIRTUAL_ENV=$PWD/../../.venv PATH=$PWD/../../.venv/bin:$PATH maturin develop --release
```

## Inputs

The default set (`common.default_inputs()`):

| Name | n | median length | total residues | Source |
|------|--:|--------------:|---------------:|--------|
| `sample36` | 36 | 377 | 13,795 | `mafft-upstream/test/sample` (protein) |
| `mtb_cds_120x1400.fa` | 120 | 1,415 | 169,859 | `crates/mafft-bin/tests/fixtures/` (DNA CDS) |
| `BB30003.tfa` | 142 | 401.5 | 57,849 | BAliBASE 3 RV30 |
| `BB30013.tfa` | 86 | 671.5 | 56,438 | BAliBASE 3 RV30 |
| `BB30001.tfa` | 116 | 435 | 54,289 | BAliBASE 3 RV30 |
| `BB30005.tfa` | 113 | 491 | 54,215 | BAliBASE 3 RV30 |
| `BB40049.tfa` | 62 | 778 | 53,444 | BAliBASE 3 RV40 |

The five BAliBASE sets are the largest of the 386 in BAliBASE 3 by total
(ungapped) residues. BAliBASE is not redistributed here; the drivers look
for it under `~/data/balibase/bb3_release` or `$BALIBASE_DIR` and skip
missing files. Any FASTA path can be passed with `--inputs`; the
determinism-only breadth sweep uses all 76 sets of BAliBASE 3 RV11
(4–14 sequences, median of per-set median lengths 259, 143,654 residues in
total).

## Running

```sh
# Peak RSS: 4 modes x --thread 1,4,8 x 3 runs, plus --memsave on the two largest inputs
python3 benchmarks/memory.py
python3 benchmarks/summarize_memory.py benchmarks/results/memory_<date>.tsv [later.tsv ...]
# Check that timing the wrapper covers every engine binary it spawns
python3 benchmarks/memory.py --verify-wrapper --inputs sample36 --modes fftnsi linsi --threads 1 4

# Determinism + speedup (FFT-NS-i, L-INS-i, G-INS-i; --thread 1,2,4,8,10; 3 interleaved runs)
python3 benchmarks/threads.py --inputs sample36 mtb_cds_120x1400.fa BB40049.tfa --max-load 3
# Determinism-only sweep over a whole BAliBASE subset, mafft-rs only
python3 benchmarks/threads.py --impls rust --threads 1 4 8 \
    --inputs ~/data/balibase/bb3_release/RV11/*.tfa --out benchmarks/results/threads_rv11_determinism.tsv
# Re-print the tables of an existing TSV
python3 benchmarks/threads.py --summarize benchmarks/results/threads_<date>.tsv

# pymafft threads= and GIL release
.venv/bin/python benchmarks/pymafft_threads.py
```

## How the numbers are taken

- **Peak RSS** is the "maximum resident set size" of `/usr/bin/time -l`
  (macOS, bytes) or `/usr/bin/time -v` (Linux, KiB, converted). For C the
  timed process is the `mafft` wrapper script; the kernel reports the
  maximum over the process and its waited-for descendants, so the figure
  is the peak of the largest engine binary of the pipeline (`disttbfast`,
  `tbfast`, `dvtditr`, …), not a sum. `memory.py --verify-wrapper`
  checks this by shimming every engine binary with its own
  `/usr/bin/time` and comparing the per-binary maximum with the wrapper's
  figure (they were equal on macOS 26 / Apple M4 for FFT-NS-i and L-INS-i
  at `--thread 1` and `--thread 4`).
- **Wall clock** is measured around the subprocess in Python
  (`time.perf_counter`), so it includes process start-up and, for C, the
  wrapper script and its temporary files — what a user of either command
  waits for. C and `mafft-rs` runs alternate (C, rs, rs, C, …) and the
  median is reported.
- **Determinism** is a SHA-256 of stdout. `threads.py` exits non-zero if
  any `mafft-rs` output differs from `mafft-rs --thread 1`.
- `--thread N` is passed to both sides at every thread count. In C,
  `--thread 1` and no `--thread` are different code paths (the `athread`
  refinement rules); `mafft-rs` follows the same split, so the
  `--thread 1` output is the byte-identity reference here.

## Results in this directory (Apple M4, 10 cores, 32 GB, macOS 26.5)

| File | What |
|------|------|
| `results/memory_2026-10-31.tsv` | full `memory.py` grid, 522 runs, `mafft-rs` at `8bff99f` |
| `results/memory_parttree_rust_2026-10-31.tsv` | `mafft-rs --parttree` re-measured after the PartTree port (#44); at `8bff99f` it panicked above 50 sequences |
| `results/memory_2026-10-31_summary.md` | per-cell medians of the two files above |
| `results/threads_2026-10-31.tsv` | `threads.py` timing + determinism: sample36 (FFT-NS-i, L-INS-i, G-INS-i), 120 × 1.4 kb CDS (FFT-NS-i), BB30005 (all three); load 2.2–5.9 |
| `results/threads_bb_largest_determinism_2026-10-31.tsv` | `mafft-rs` only, five largest BAliBASE sets, FFT-NS-i / L-INS-i, `--thread 1 4 8`, 2 runs (timings not meaningful: ran alongside `memory.py`) |
| `results/threads_rv11_determinism_2026-10-31.tsv` | `mafft-rs` only, all 76 RV11 sets, FFT-NS-i / L-INS-i / G-INS-i, `--thread 1 4 8`, 3 runs (timings not meaningful) |
| `results/pymafft_threads_2026-10-31.tsv` | `pymafft_threads.py`: sample36 and BB30005, plus the GIL demonstration |

Interpretation and the headline tables: `PROFILING.md` ("Memory
footprint", "Multi-threading") and the README "Performance" section.

## Caveats

- Wall-clock results depend on machine load. The drivers record the
  1-minute load average with every row (`threads.py`, `pymafft_threads.py`)
  and `threads.py --max-load X` waits for the load to drop before each
  batch. Peak RSS and output hashes do not depend on load.
- macOS reports `ru_maxrss` in bytes and Linux in KiB; `common.parse_rss`
  converts both to bytes. The `peak memory footprint` line that macOS
  `time -l` also prints covers only the direct child and is not used.
- C's `--thread N` for N >= 2 is nondeterministic, so its timing at
  those thread counts is for whichever alignment that run produced.
