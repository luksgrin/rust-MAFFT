//! FFT-based homology detection for MAFFT.
//!
//! Uses `rustfft` with `num_complex::Complex64` (replacing the C hand-rolled
//! Cooley-Tukey FFT on `Fukusosuu` arrays).
//!
//! The three main operations are:
//! 1. **Cross-correlation** — FFT-based correlation between two sequence
//!    profiles to find optimal lag (shift) positions.
//! 2. **Segment detection** — sliding-window scoring to identify alignable
//!    regions above a threshold.
//! 3. **Candidate selection** — picking the top-k lag positions from the
//!    correlation profile.

mod correlation;
mod segments;
mod candidates;

pub use correlation::{cross_correlate, inner_product};
pub use segments::{alignable_segments, AlignableSegment, SegmentParams};
pub use candidates::get_top_candidates;
