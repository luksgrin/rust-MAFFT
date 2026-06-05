# mafft-io

> FASTA, Clustal, PHYLIP, hat2, and localhom I/O for [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT).

> **Acknowledgment.** rust-MAFFT is a port of [MAFFT](https://mafft.cbrc.jp/alignment/software/)
> by **Kazutaka Katoh** and colleagues (CBRC). The science is theirs.
> If you use this crate in published work, please cite
> [Katoh & Standley 2013](https://doi.org/10.1093/molbev/mst010);
> full guidance at [the project citation page](https://luksgrin.github.io/rust-MAFFT/citation/).

Reads sequences and their associated constraint / distance data from
disk (or any `BufRead`), and writes alignment results back out. Reuses
[`noodles-fasta`](https://crates.io/crates/noodles-fasta) under the
hood for the FASTA parser. Handles MAFFT's case-preservation and
gap-character normalisation conventions so the engine sees byte-for-
byte the same input as the C reference.

## Install

```sh
cargo add mafft-io
```

Most callers want the top-level [`mafft`](https://crates.io/crates/mafft)
crate, which re-exports the public surface here.

## Examples

```rust
use mafft_io::{read_fasta, read_fasta_from_reader};

// From a path
let input = read_fasta("input.fasta")?;

// From any BufRead (stdin, a Vec<u8>, an HTTP body, ...)
let cursor = std::io::Cursor::new(b">a\nACDEF\n>b\nACDF\n");
let input  = read_fasta_from_reader(std::io::BufReader::new(cursor))?;
```

## Documentation

API reference: [docs.rs/mafft-io](https://docs.rs/mafft-io)

Project documentation: [luksgrin.github.io/rust-MAFFT](https://luksgrin.github.io/rust-MAFFT/)

## License

MIT for the Rust port; the upstream MAFFT I/O conventions this crate
mirrors are BSD-3-Clause. See `LICENSE-MIT` and `LICENSE-BSD` in the
workspace root.
