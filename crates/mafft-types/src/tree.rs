// Tree types are defined in the `mafft-tree` crate (`Topology`, `JoinStep`)
// rather than here, because they are tightly coupled to the tree-building
// algorithms. This module is intentionally empty.
//
// The C `Treedep` struct is represented by `mafft-tree::Topology` which
// tracks child relationships and branch lengths in its `JoinStep` entries.
