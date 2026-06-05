# mafft-rs

> The MAFFT multiple sequence alignment CLI, in pure Rust. Drop-in replacement for upstream `mafft`.

Byte-identical output to C MAFFT 7.526 across the full BAliBASE 3
fixture set (1930/1930). Same flag surface (`--localpair`,
`--globalpair`, `--genafpair`, `--maxiterate`, `--add`,
`--adjustdirection`, `--parttree`, …), same FASTA output format, same
exit codes. The only visible difference is the binary name —
`mafft-rs` instead of `mafft`.

## Install

```sh
# Via cargo (requires Rust toolchain)
cargo install mafft-rs

# Via pip (no Rust toolchain required; the wheel bundles the binary)
pip install pymafft

# Or grab a pre-built binary from
# https://github.com/luksgrin/rust-MAFFT/releases
```

## Usage

```sh
# Default (FFT-NS-2)
mafft-rs input.fasta > aligned.fasta

# Most accurate for < 200 sequences
mafft-rs --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Suppress progress
mafft-rs --quiet input.fasta > aligned.fasta
```

## Documentation

CLI reference: [luksgrin.github.io/rust-MAFFT/cli-reference/](https://luksgrin.github.io/rust-MAFFT/cli-reference/)

Project documentation: [luksgrin.github.io/rust-MAFFT](https://luksgrin.github.io/rust-MAFFT/)

For the programmatic library, see
[`mafft`](https://crates.io/crates/mafft) or
[`mafft-core`](https://crates.io/crates/mafft-core).

## License

MIT for the Rust port; the algorithms are derived from upstream MAFFT
under BSD-3-Clause. See `LICENSE-MIT` and `LICENSE-BSD` in the
workspace root.
