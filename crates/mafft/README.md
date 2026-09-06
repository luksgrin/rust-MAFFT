# mafft

MAFFT multiple sequence alignment — high-level Rust API.

This is the ergonomic entry point to the [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT)
workspace. It re-exports the alignment engine, sequence types, and I/O
helpers so most callers can `use mafft::*` and skip the sub-crate
imports.

> **Acknowledgment.** rust-MAFFT is a port of [MAFFT](https://mafft.cbrc.jp/alignment/software/)
> by **Kazutaka Katoh** and colleagues (CBRC). The science is theirs.
> If you use this crate in published work, please cite
> [Katoh & Standley 2013](https://doi.org/10.1093/molbev/mst010);
> full guidance at [the project citation page](https://luksgrin.github.io/rust-MAFFT/citation/).

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

## The command line's flag layer (`cli` feature)

`--auto` chooses the mode from the input; `--adjustdirection` flips strands
before aligning; `--nuc` / `--amino` force the type and the case fold. To
get exactly that behaviour in-process instead of re-deriving it, enable the
`cli` feature and use `mafft::cli`, which re-exports the
[`mafft-rs`](https://crates.io/crates/mafft-rs) library:

```toml
mafft = { version = "0.1", features = ["cli"] }
```

```rust,no_run
use mafft::cli::{run_from_seqs, SilentProgress};
use mafft::{read_fasta};

let input = read_fasta("input.fasta").unwrap();
let msa = run_from_seqs(["mafft", "--auto", "--nuc"], &input, &SilentProgress).unwrap();
```

## When to depend on a sub-crate instead

Use the underlying crates directly only when you need finer control:

- `mafft-core` — engine + modes (skip the I/O layer)
- `mafft-align`, `mafft-tree`, `mafft-scoring`, `mafft-fft` — internals

## Related

- [`mafft-rs`](https://crates.io/crates/mafft-rs) — standalone CLI: `cargo install mafft-rs`
- [`pymafft`](https://pypi.org/project/pymafft/) — Python bindings: `pip install pymafft`

## License

MIT AND BSD-3-Clause (the upstream MAFFT C source is BSD-licensed).
