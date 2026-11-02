# High-throughput workload: full run (issue #23)

This is the hand-over document for the *full* timed run of the
high-throughput workload. The corpus, drivers and a 50-cluster smoke run
were merged in #42; what remains is the measurement itself, which takes
about a day per replicate and therefore has to be scheduled on a machine
that nobody else is using. Everything below is self-contained: a person
who has never seen this repository should be able to produce the results
TSV and the summary table from it.

## What is being measured, and why

Pipelines such as panaroo, Roary or OrthoFinder call MAFFT once per gene
cluster, thousands of times per run. Each call pays for a process spawn,
FASTA serialisation of the input, and parsing of the output. The claim
under test is that calling `pymafft.align()` in-process removes that
overhead, and the question is how much that is worth at the scale of a
real pipeline, on one core and on several.

The same 2,695 DNA clusters (and optionally 1,500 protein orthogroups)
are aligned seven ways with identical flags
(`--auto --adjustdirection --thread 1 --nuc`, the panaroo command line):

| `--impl`     | what it is                                                        | what it isolates |
|--------------|-------------------------------------------------------------------|------------------|
| `snake-c`    | Snakemake, one job per cluster, C MAFFT (in-tree build)           | the baseline a pipeline author has today |
| `snake-rust` | same Snakefile, `mafft-rs`                                        | engine speed, same scheduler |
| `loop-c`     | plain Python thread pool spawning the C wrapper per cluster       | spawn + FASTA cost without Snakemake's scheduler |
| `loop-rust`  | same loop, `mafft-rs`                                             | as above, Rust engine |
| `py-serial`  | one process, `pymafft.align()` per cluster                        | in-process, no parallelism |
| `py-process` | `pymafft.align()` in a `ProcessPoolExecutor`                      | in-process, N workers |
| `py-thread`  | `pymafft.align()` in a `ThreadPoolExecutor`                       | in-process, N threads; equal to `py-process` only if the GIL is really released |

Each cell records total wall-clock, user+system CPU, peak RSS of the whole
process tree (sampled), the number of clusters aligned, and whether every
output is byte-identical to the `snake-c` output. Workers are 1, 4 and 8.

## Hardware and conditions

- Any x86-64 or arm64 machine with ≥ 8 physical cores and ≥ 16 GB RAM.
  The published smoke run and all other benchmarks of this repository are
  from an Apple M4 (10 cores, 32 GB); a Linux x86-64 run is explicitly
  welcome — it is the other reference platform
  (see `docs/architecture/byte-identity.md`).
- **The machine must be otherwise idle for the whole run.** Load from
  anything else (another benchmark, a browser, a backup) invalidates the
  wall-clock numbers. Check `uptime` before starting; the load average
  should be below 1.
- **Prevent sleep.** On macOS prefix every long command with
  `caffeinate -i`; on Linux make sure no idle suspend is configured.
- Disk: ~1 GB for the corpus and scratch output (`target/`).
- Budget: C needs ~2 s CPU per DNA cluster and ~6 s per orthogroup
  (`--auto` selects L-INS-i for inputs of this size). One single-worker C
  run over the DNA set is ~1.6 h; the default DNA grid is 17 runs
  (7 implementations × workers, minus `py-serial` at > 1 worker), roughly
  **one day for DNA with one replicate**, plus ~1.5 days if protein is
  included. The run is resumable (cached per cell), so it can be split
  over several nights.

## One-time setup

All commands run from a checkout of `luksgrin/rust-MAFFT` at the commit
you want to measure (record it: `git rev-parse HEAD`). Use a fresh clone,
not a worktree shared with other work.

```sh
git clone --recurse-submodules https://github.com/luksgrin/rust-MAFFT.git
cd rust-MAFFT
git rev-parse HEAD                     # record this

# 1. Rust toolchain (https://rustup.rs) and the two binaries
cargo build --release -p mafft-rs      # target/release/mafft-rs

# 2. The C reference, built in-tree from the pinned submodule.
#    Do NOT use a downloaded/bioconda mafft: it is a different build and
#    byte-identity is defined against the in-tree one.
make -C mafft-upstream/core            # -> mafft-upstream/binaries, scripts/mafft
export MAFFT_BINARIES="$PWD/mafft-upstream/binaries"
mafft-upstream/scripts/mafft --version # must print v7.526

# 3. Python side (3.9-3.14; 3.12 is what the smoke run used)
python3.12 -m venv .venv
.venv/bin/pip install --upgrade pip maturin biopython snakemake
( cd crates/pymafft && ../../.venv/bin/maturin develop --release )
.venv/bin/python -c "import pymafft, snakemake; print(pymafft.__version__)"

# 4. The corpus: re-downloads the exact committed accessions and verifies
#    checksums. Needs network; ~5-10 min. MMseqs2 is NOT needed (only for
#    --resample, which you must not use).
.venv/bin/python benchmarks/htworkload_corpus.py build
.venv/bin/python benchmarks/htworkload_corpus.py verify   # must say OK
.venv/bin/python benchmarks/htworkload_corpus.py table    # composition, for the report
```

If `build` fails on a transient HTTP error, run it again; it resumes.

## Sanity check before the long run (10 minutes)

```sh
uptime                                  # load < 1
caffeinate -i .venv/bin/python benchmarks/htworkload.py run --smoke --molecule dna
.venv/bin/python benchmarks/htworkload.py summarize benchmarks/results/htworkload_<os>-<arch>_<date>.tsv
```

Expect: all 50 clusters aligned on every path and `identical = yes` for
every Rust/pymafft cell. Cluster `dna_01715` diverged in the committed
smoke TSV (issue #43); that FFT-path divergence is fixed and the cluster
is now a byte-identity fixture. If anything is not identical, stop and
report it before measuring.

## The full run

```sh
# DNA, the primary result. Resumable: re-running fills in missing cells.
caffeinate -i .venv/bin/python benchmarks/htworkload.py run --molecule dna

# Optional second night: protein orthogroups
caffeinate -i .venv/bin/python benchmarks/htworkload.py run --molecule protein

# Optional third+ nights: replicates (adds rep 2 and 3 of every cell)
caffeinate -i .venv/bin/python benchmarks/htworkload.py run --molecule dna --reps 3
```

Notes

- `run --list` prints the cells without running anything; `--impl` and
  `--workers` restrict the grid if a night is short, e.g.
  `--impl snake-c py-process --workers 1 8`.
- Do not use `--limit` or `--smoke` for the real run; those write to a
  separate cache and are for pipeline checks only.
- The results TSV is `benchmarks/results/htworkload_<os>-<arch>_<date>.tsv`;
  per-cell caches are under `target/htworkload/<host>/`. If the machine is
  interrupted, just rerun the same command.
- `--keep` keeps every aligned FASTA (several hundred MB); not needed.
- The driver records the 1-minute load average at the start of every
  cell (`load1` column of the TSV). After the run, check that column:
  any cell started under a load above ~1.5 is suspect — delete its cache
  directory (`target/htworkload/<host>/<molecule>/<impl>_w<workers>_rep<n>`)
  and rerun the same command to redo just that cell.

## Reporting

```sh
.venv/bin/python benchmarks/htworkload.py summarize benchmarks/results/htworkload_<os>-<arch>_<date>.tsv
```

prints, per implementation × workers: total wall-clock, CPU, peak RSS,
per-alignment time, speed-up relative to `snake-c` at the same worker
count, and the byte-identity count. Deliver:

1. the results TSV, committed under `benchmarks/results/`;
2. the summary table and the `htworkload_corpus.py table` output pasted
   into `PROFILING.md` under "High-throughput workload", with the host
   (CPU model, core count, RAM, OS), the commit measured, the date, and
   the load average observed;
3. a short note of any cell that is not byte-identical, with the
   cluster id, as a comment on issue #23.

Open a pull request against `main` with `Closes #23`. That PR touches only
`benchmarks/` and `PROFILING.md`, so no CI runs; that is expected.

## Interpreting the result

The quantity of interest is the ratio `snake-c` / `py-process` (and
`loop-c` / `py-process`) at each worker count: that is what a pipeline
gains by replacing the subprocess with the in-process call. The
`snake-rust` vs `snake-c` ratio separates engine speed from call overhead.
`py-thread` ≈ `py-process` confirms that the bindings release the GIL. If
`py-process` at 8 workers is not close to 8× `py-serial`, say so; the
clusters are small and process-pool overhead is itself part of the result.
