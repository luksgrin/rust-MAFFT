# mafft-tree

> Distance computation and guide tree construction for [rust-MAFFT](https://github.com/luksgrin/rust-MAFFT).

k-mer distances, full pairwise distance matrices, NJ and UPGMA tree
construction, PartTree partitioning for large datasets (10K+
sequences), and the memsavetree variants (`memsavetree`,
`compacttree`, youngest/average/minimum/mixed linkage). All
byte-identical to the C reference: the same merge order, the same
tie-breaks, the same per-node weight propagation.

## Install

```sh
cargo add mafft-tree
```

Most callers want the top-level [`mafft`](https://crates.io/crates/mafft)
crate. Depend on `mafft-tree` directly when you need guide trees for
something other than progressive MSA — e.g. building a quick rough
phylogeny from a k-mer distance matrix.

## Documentation

API reference: [docs.rs/mafft-tree](https://docs.rs/mafft-tree)

Project documentation: [luksgrin.github.io/rust-MAFFT](https://luksgrin.github.io/rust-MAFFT/)

## License

MIT for the Rust port; the tree-building algorithms are derived from
upstream MAFFT under BSD-3-Clause. See `LICENSE-MIT` and `LICENSE-BSD`
in the workspace root.
