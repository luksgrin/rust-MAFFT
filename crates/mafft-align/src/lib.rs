//! Pairwise and multi-sequence alignment algorithms for MAFFT.
//!
//! Phase 4a: Single-sequence pairwise alignment (Galign11, Lalign11, genalign11).
//! Phase 5: Profile alignment (MSalignmm), FFT-accelerated alignment (Falign),
//!          and local homology constraint building (pairlocalalign).

mod dp;
mod global;
mod local;
mod genaffine;
mod profile;
mod fft_align;
mod constraints;

pub use dp::{Alignment, AlignOp, GapModel};
pub use global::global_align;
pub use local::{local_align, LocalAlignment};
pub use genaffine::{genaffine_local_align, GenAffineGapModel};
pub use profile::{Profile, profile_align};
pub use fft_align::{fft_profile_align, FftAlignParams, Anchor};
pub use constraints::build_local_homology_table;
