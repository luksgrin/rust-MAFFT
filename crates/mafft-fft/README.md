# mafft-fft

> Bit-for-bit-compatible Cooley-Tukey FFT for [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT).

> **Acknowledgment.** rust-MAFFT is a port of [MAFFT](https://mafft.cbrc.jp/alignment/software/)
> by **Kazutaka Katoh** and colleagues (CBRC). The science is theirs.
> If you use this crate in published work, please cite the foundational
> FFT-based MSA paper, [Katoh et al. 2002](https://doi.org/10.1093/nar/gkf436),
> together with [Katoh & Standley 2013](https://doi.org/10.1093/molbev/mst010);
> full guidance at [the project citation page](https://luksgrin.github.io/rust-MAFFT/citation/).

A hand port of MAFFT's `core/fft.c`. Drives the FFT-anchored homology
detection used by FFT-NS-2 / FFT-NS-i and the Falign-routed iterative
modes. Matches the C reference rounding exactly — off-the-shelf FFT
crates (rustfft, realfft, FFTW) all vary in butterfly grouping order,
and on flat-landscape correlation matrices those 1-ULP differences
flip anchor selection and ultimately the alignment.

## Install

```sh
cargo add mafft-fft
```

Most callers should depend on [`mafft`](https://crates.io/crates/mafft)
or [`mafft-core`](https://crates.io/crates/mafft-core) and let them
pull in `mafft-fft` transitively. Direct use is only useful if you're
building an alignment-adjacent tool that wants the same FFT primitives
without the rest of the engine.

## Documentation

API reference: [docs.rs/mafft-fft](https://docs.rs/mafft-fft)

Project documentation: [luksgrin.github.io/rust-MAFFT](https://luksgrin.github.io/rust-MAFFT/) ·
[Byte-identity design note](https://luksgrin.github.io/rust-MAFFT/architecture/byte-identity/)

## License

MIT for the Rust port; the FFT algorithm matches upstream MAFFT under
BSD-3-Clause. See `LICENSE-MIT` and `LICENSE-BSD` in the workspace
root.
