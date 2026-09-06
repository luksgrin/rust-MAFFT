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

    Verify:

    ```rust
    use mafft::{MafftEngine, AlignmentMode, SequenceSet, Sequence};

    let input = SequenceSet {
        sequences: vec![
            Sequence { name: "a".into(), data: b"ACDEFGHIK".to_vec() },
            Sequence { name: "b".into(), data: b"ACDEFHIK".to_vec() },
        ],
        seq_type: mafft::SeqType::Protein,
    };
    let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
    assert_eq!(msa.sequences.len(), 2);
    ```

=== "Python package"

    ```sh
    pip install pymafft
    ```

    Wheels are published for Linux (x86_64, aarch64), macOS (Intel,
    Apple Silicon), and Windows (x86_64) across Python 3.9 – 3.13. The
    install also drops a `mafft-rs` console script on your `$PATH` —
    the same native binary, just reachable from any pip environment.

    Verify:

    ```python
    import pymafft
    result = pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])
    print(result.to_fasta())
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
    mafft-rs sample.fasta > aligned.fasta
    ```

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
cargo test --workspace --exclude pymafft --release   # full test suite (537 tests + 7 ignored)
```

The `mafft-upstream/` submodule pins the reference C MAFFT version
that the FFI cross-validation harness builds against. It's only
required for development — release builds compile zero C code.
