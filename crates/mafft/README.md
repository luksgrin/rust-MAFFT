# mafft

MAFFT multiple sequence alignment — high-level Rust API.

This is the ergonomic entry point to the [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT)
workspace. It re-exports the alignment engine, sequence types, and I/O
helpers so most callers can `use mafft::*` and skip the sub-crate
imports.

## Install

```sh
cargo add mafft
```

## Usage

```rust,no_run
use mafft::{MafftEngine, AlignmentMode, read_fasta};

let input = read_fasta("input.fasta").unwrap();
let engine = MafftEngine::new(AlignmentMode::FftNs2);
let msa = engine.align(&input);
```

## Available strategies

| Mode | Use case |
|------|----------|
| `AlignmentMode::FftNs2` | Default — fast progressive alignment |
| `AlignmentMode::FftNsi` | Iterative refinement |
| `AlignmentMode::LInsi` | Most accurate for < 200 sequences |
| `AlignmentMode::GInsi` | Globally alignable sequences |
| `AlignmentMode::EInsi` | Sequences with large gaps |

## When to depend on a sub-crate instead

Use the underlying crates directly only when you need finer control:

- `mafft-core` — engine + modes (skip the I/O layer)
- `mafft-align`, `mafft-tree`, `mafft-scoring`, `mafft-fft` — internals

## Related

- [`mafft-rs`](https://crates.io/crates/mafft-rs) — standalone CLI: `cargo install mafft-rs`
- [`pymafft`](https://pypi.org/project/pymafft/) — Python bindings: `pip install pymafft`

## License

MIT AND BSD-3-Clause (the upstream MAFFT C source is BSD-licensed).
