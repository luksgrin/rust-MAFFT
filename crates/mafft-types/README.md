# mafft-types

> Pure-Rust types shared across the [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT) workspace.

`Sequence`, `SequenceSet`, `SeqType`, `ScoringModel`, gap-segment
descriptors, `HomologyRegion` for local-homology constraint lists,
plus a `Complex64` alias used by the FFT layer. No engine logic, no
I/O — just the data structures every other rust-MAFFT crate agrees
on.

## Install

```sh
cargo add mafft-types
```

If you're writing an application, you almost certainly want the
top-level [`mafft`](https://crates.io/crates/mafft) crate instead,
which re-exports the public types from here. Use `mafft-types`
directly only when you need a minimal-dependency way to construct
inputs or interpret outputs for code that consumes the engine via a
different crate.

## Documentation

API reference: [docs.rs/mafft-types](https://docs.rs/mafft-types)

Project documentation: [luksgrin.github.io/rust-MAFFT](https://luksgrin.github.io/rust-MAFFT/)

## License

MIT for the Rust port; the upstream MAFFT algorithmic constructs this
crate's types model are BSD-3-Clause. See `LICENSE-MIT` and
`LICENSE-BSD` in the workspace root.
