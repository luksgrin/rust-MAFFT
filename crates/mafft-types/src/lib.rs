//! Pure Rust equivalents of MAFFT's core C data structures.
//!
//! These types use idiomatic Rust (Vec, Option, enums) instead of raw pointers
//! and linked lists. Conversion traits to/from the C representations in
//! `mafft_sys` are provided for the FFI boundary.

mod local_hom;
mod tree; // intentionally empty — tree types live in mafft-tree
mod segment;
mod complex;
mod seq;
mod scoring;

pub use local_hom::*;
pub use segment::*;
pub use complex::*;
pub use seq::*;
pub use scoring::*;
