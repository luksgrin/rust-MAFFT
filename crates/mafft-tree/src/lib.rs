//! Distance computation and guide tree construction for MAFFT.
//!
//! Ports the C `nj.c` (neighbor-joining), UPGMA variants from `mltaln9.c`,
//! k-tuple distance from `tddis.c`, and sequence weighting from tree topology.

mod distance;
mod nj;
mod upgma;
mod musclesupg;
mod addonetip;
pub mod parttree_dist;
pub mod parttree_pivot;
pub mod parttree_split;
mod topology;
mod weighting;
pub mod bsd_qsort;
pub mod newick;
pub mod treein;
pub mod memsavetree;

pub use distance::{DistanceMatrix, pairwise_identity_distance, ktuple_distance, scoring_matrix_distance};
pub use nj::neighbor_joining;
pub use upgma::{upgma, upgma_int};
pub use musclesupg::{musclesupg, ClusterMethod};
pub use addonetip::{addonetip, compute_distfromtip, AddResult};
pub use topology::{Topology, JoinStep};
pub use weighting::{sequence_weights, BranchWeights};
pub use newick::topology_to_newick;
pub use treein::{parse_mafft_tree, parse_mafft_tree_str};
