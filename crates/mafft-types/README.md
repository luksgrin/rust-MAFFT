# mafft-types

> Pure-Rust types shared across the [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT) workspace.

> **Acknowledgment.** rust-MAFFT is a port of [MAFFT](https://mafft.cbrc.jp/alignment/software/)
> by **Kazutaka Katoh** and colleagues (CBRC). The science is theirs.
> If you use this crate in published work, please cite
> [Katoh & Standley 2013](https://doi.org/10.1093/molbev/mst010);
> full guidance at [the project citation page](https://luksgrin.github.io/rust-MAFFT/citation/).

`Sequence`, `SequenceSet`, `SeqType`, `ScoringModel`, gap-segment
descriptors, `HomologyRegion` for local-homology constraint lists,
plus a `Complex64` alias used by the FFT layer. No engine logic, no
I/O — just the data structures every other rust-MAFFT crate agrees
on, and the floating-point contraction policy `mafft_types::fp`
(`CONTRACTS_FMA`, `fmadd`; cargo features `fp-contract-fma` /
`fp-contract-none`) that decides whether `a*b + c` rounds once or
twice, mirroring the reference C build of the target platform.

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
