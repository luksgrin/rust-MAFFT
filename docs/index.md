# rust-MAFFT

**Pure-Rust port of [MAFFT 7.526](https://mafft.cbrc.jp/alignment/software/), byte-identical to the C reference.**

!!! quote "Built on MAFFT"

    rust-MAFFT is a port — every algorithmic decision in this codebase
    traces back to work by **[Kazutaka Katoh](https://mafft.cbrc.jp/alignment/software/)**
    and colleagues at CBRC. The FFT-anchored alignment algorithm, the
    iterative-refinement strategies, the scoring matrices, the PartTree
    heuristic — all of it is theirs. This project is engineering on
    top of their science. **If you use it in published work, please
    cite the original MAFFT paper.** See the
    [Citation](citation.md) page for full details and BibTeX.

## What it is

A complete rewrite of the MAFFT multiple sequence alignment toolchain in
Rust — no `unsafe`, no C at runtime, no global state. Output matches C
MAFFT 7.526 byte-for-byte across BAliBASE 3 (1930/1930 fixtures) for
every supported mode: FFT-NS-2, FFT-NS-i, G-INS-i, L-INS-i, E-INS-i,
PartTree, `--add`, `--adjustdirection`, and more.

## Three ways to use it

<div class="grid cards" markdown>

-   :material-language-rust: **As a Rust library**

    ---

    ```sh
    cargo add mafft
    ```

    Programmatic access from any Rust project. See the
    [Rust API on docs.rs](https://docs.rs/mafft).

-   :material-language-python: **As a Python package**

    ---

    ```sh
    pip install pymafft
    ```

    Native bindings via PyO3. Drop-in Biopython interop.
    See the [Python API](python-api/index.md).

-   :material-console: **As a CLI binary**

    ---

    ```sh
    # via cargo
    cargo install mafft-rs

    # or via pip (also bundles the binary)
    pip install pymafft

    # or pre-built from GitHub releases
    ```

    See the [CLI Reference](cli-reference.md).

</div>

## Why byte-identity matters

MAFFT is the *de facto* MSA reference for bioinformatics workflows. A
Rust port that drifts even one tie-break from the C output forces every
downstream pipeline to choose: re-validate the rewrite, or maintain a
parallel C dependency. Byte-identity removes that choice — the Rust
binary is a drop-in replacement.

The [Byte-identity](architecture/byte-identity.md) page documents the
design decisions, FFT tie-break ports, FP-order constraints, and
cross-validation infrastructure (99 FFI test functions against the C
implementation, compiled in-tree on the same platform).

## Status

- **All standard modes**: byte-identical across BAliBASE 3.
- **Stub-parity flags** (`--pdbidlist`, `--pdbfilelist`): match C MAFFT's
  "temporarily unavailable, 2018/Dec." exit verbatim — these were
  disabled upstream.
- **DNA**: byte-identical on nucleotide input too (pair-phase gap scale,
  `dndpre` offset, lowercase output, `--nuc` / `--amino`), including C's
  distinct `--thread 1` refinement path.
- **In-memory API**: `mafft_rs::run_from_seqs` / `Mafft::run_seqs` and
  `pymafft.align(...)` share the CLI's flag layer, so `--auto`,
  `--adjustdirection`, `--nuc` behave identically from Rust, Python and
  the shell.
- **Out of scope / residuals**: see the [architecture index](architecture/index.md).

## Quick links

- :material-source-branch: [GitHub](https://github.com/luksgrin/rust-MAFFT)
- :material-package: [crates.io: `mafft`](https://crates.io/crates/mafft) · [`mafft-rs`](https://crates.io/crates/mafft-rs) · [`mafft-core`](https://crates.io/crates/mafft-core)
- :material-language-python: [PyPI: `pymafft`](https://pypi.org/project/pymafft/)
- :material-book-open-page-variant: [docs.rs: `mafft`](https://docs.rs/mafft)
- :material-format-quote-close: [**How to cite**](citation.md)
