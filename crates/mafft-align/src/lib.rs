//! Pairwise and multi-sequence alignment algorithms for MAFFT.
//!
//! Phase 4a: Single-sequence pairwise alignment (Galign11, Lalign11, genalign11).
//! Phase 5: Profile alignment (MSalignmm), FFT-accelerated alignment (Falign),
//!          local homology constraint building (pairlocalalign),
//!          and constrained alignment (Falign_localhom, partSalignmm).

mod dp;
mod global;
mod local;
mod genaffine;
mod profile;
mod fft_align;
mod constraints;
mod constrained_align;

pub use dp::{Alignment, AlignOp, GapModel};
pub use global::global_align;
pub use local::{local_align, LocalAlignment};
pub use genaffine::{genaffine_local_align, GenAffineGapModel};
pub use profile::{
    Profile, profile_align, profile_align_imp, profile_align_imp_with_tiebreak,
    profile_align_imp_with_boundary, BoundaryFreqs, pairwise_align11, align_with_anchors,
};
pub use fft_align::{fft_profile_align, find_fft_anchors, FftAlignParams, Anchor};
pub use constraints::{
    build_local_homology_table, build_homology_table,
    build_homology_table_with_unalign, build_imp_matrix,
    recompute_importance, PairAligner, FASTATHRESHOLD_DEFAULT,
};
pub use constrained_align::{
    constrained_profile_align, partial_profile_align, ConstrainedAlignParams,
};
