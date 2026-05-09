//! FFT-based homology detection for MAFFT.
//!
//! Uses `rustfft` with `num_complex::Complex64` (replacing the C hand-rolled
//! Cooley-Tukey FFT on `Fukusosuu` arrays).

mod correlation;
mod segments;
mod candidates;
mod vectorize;
mod block_align;
pub mod fft_c_compat;

pub use correlation::{cross_correlate, inner_product};
pub use segments::{alignable_segments, AlignableSegment, SegmentParams};
pub use candidates::get_top_candidates;
pub use vectorize::{
    sequences_to_channels, sequences_to_property_channels,
    multichannel_correlate, PROTEIN_CHANNELS, DNA_CHANNELS, PROPERTY_CHANNELS,
};
pub use block_align::{block_align, block_align3};
