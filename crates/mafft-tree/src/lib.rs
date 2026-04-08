//! Distance computation and guide tree construction for MAFFT.
//!
//! Ports the C `nj.c` (neighbor-joining), UPGMA variants from `mltaln9.c`,
//! k-tuple distance from `tddis.c`, and sequence weighting from tree topology.

mod distance;
mod nj;
mod upgma;
mod musclesupg;
mod topology;
mod weighting;

pub use distance::{DistanceMatrix, pairwise_identity_distance, ktuple_distance};
pub use nj::neighbor_joining;
pub use upgma::{upgma, upgma_int};
pub use musclesupg::{musclesupg, ClusterMethod};
pub use topology::{Topology, JoinStep};
pub use weighting::sequence_weights;
