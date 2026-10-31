# Install

rust-MAFFT ships through three independent channels. Pick whichever
fits your workflow — they all produce byte-identical alignments.

=== "Rust library"

    Add the top-level `mafft` crate to your project:

    ```sh
    cargo add mafft
    ```

    `mafft` is a thin re-export of `mafft-core`, `mafft-types`, and
    `mafft-io` — most callers should depend on it and stay stable
    across workspace reorganisations. For finer control, depend on the
    sub-crates directly.

    Verify with this `src/main.rs` and `cargo run`:

    ```rust
    use mafft::{AlignmentMode, MafftEngine, SeqType, Sequence, SequenceSet};

    fn main() {
        let input = SequenceSet {
            sequences: vec![
                Sequence { name: "a".into(), data: b"ACDEFGHIK".to_vec() },
                Sequence { name: "b".into(), data: b"ACDEFHIK".to_vec() },
                Sequence { name: "c".into(), data: b"ACDHIK".to_vec() },
            ],
            seq_type: SeqType::Protein,
        };
        let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
        for (name, row) in msa.names.iter().zip(msa.sequences.iter()) {
            println!(">{name}\n{}", std::str::from_utf8(row).unwrap());
        }
    }
    ```

=== "Python package"

    ```sh
    pip install pymafft
    ```

    Wheels are published for Linux (x86_64, aarch64), macOS (Intel,
    Apple Silicon), and Windows (x86_64). Each is a stable-ABI (`abi3`)
    wheel, so the same file installs on CPython 3.9 and every later
    version, 3.14 included. The
    install also drops a `mafft-rs` console script on your `$PATH` —
    the same native binary, just reachable from any pip environment.

    Verify:

    ```sh
    python -c "import pymafft; print(pymafft.align(['ACDEFGHIK', 'ACDEFHIK', 'ACDHIK']).to_fasta())"
    mafft-rs --version
    ```

=== "CLI binary"

    Three ways to get the standalone `mafft-rs` binary, in order of
    convenience:

    ```sh
    # Option A: cargo (requires Rust toolchain)
    cargo install mafft-rs

    # Option B: pip (no Rust toolchain needed; also installs the
    # Python bindings as a bonus)
    pip install pymafft

    # Option C: pre-built binary from the GitHub release page
    # https://github.com/luksgrin/rust-MAFFT/releases
    ```

    Verify:

    ```sh
    mafft-rs --version
    printf '>a\nACDEFGHIK\n>b\nACDEFHIK\n>c\nACDHIK\n' | mafft-rs --quiet
    ```

## Verify your installation

The checks under each tab above take seconds and need neither a C
compiler nor the `mafft-upstream` submodule. Each prints the same
three-row alignment:

```
>a
ACDEFGHIK
>b
ACDEF-HIK
>c
ACD---HIK
```

(`pymafft.align` on bare strings names the rows `seq_1`, `seq_2`,
`seq_3`.) That is all an install needs. Re-checking byte-identity
against C MAFFT is a contributor task with its own prerequisites; see
[Reproduce the cross-validation suite](../contributing/index.md#reproduce-the-cross-validation-suite).

## Platform support

| Channel | Linux x86_64 | Linux aarch64 | macOS Intel | macOS Apple Silicon | Windows x86_64 |
|--|--|--|--|--|--|
| Rust crate | ✓ | ✓ | ✓ | ✓ | ✓ |
| Python wheel | ✓ | ✓ | ✓ | ✓ | ✓ |
| Pre-built binary | ✓ | ✓ | ✓ | ✓ | ✓ |

## From source

```sh
git clone --recurse-submodules https://github.com/luksgrin/rust-MAFFT
cd rust-MAFFT
cargo build --release -p mafft-rs   # CLI
cargo test -p mafft-rs --lib        # quick check: unit tests, no C compiler needed
```

The full test suite (`cargo test --workspace --exclude pymafft --release`,
537 tests + 7 ignored) compiles the C reference and takes minutes; see
[Reproduce the cross-validation suite](../contributing/index.md#reproduce-the-cross-validation-suite).

The `mafft-upstream/` submodule pins the reference C MAFFT version
that the FFI cross-validation harness builds against. It's only
required for development — release builds compile zero C code.

The Python bindings build into a virtual environment with the standard
library's `venv` and `pip`:

```sh
cd crates/pymafft
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install maturin
maturin develop --release
```

See [Contributing](../contributing/index.md#python-bindings) for the test
dependencies and the bundled-CLI step.
