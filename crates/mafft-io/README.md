# mafft-io

> FASTA, Clustal, PHYLIP, hat2, and localhom I/O for [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT).

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
