//! Pairwise sequence alignment algorithms for MAFFT.
//!
//! Ports the C `Galign11.c` (global/Needleman-Wunsch), `Lalign11.c`
//! (local/Smith-Waterman), and `genalign11.c` (generalized affine gap).
//!
//! All algorithms use affine gap penalties with owned scoring contexts
//! instead of C globals.

mod dp;
mod global;
mod local;
mod genaffine;

pub use dp::{Alignment, AlignOp, GapModel};
pub use global::global_align;
pub use local::{local_align, LocalAlignment};
pub use genaffine::genaffine_local_align;
