//! Core MAFFT alignment engine.
//!
//! Orchestrates the full alignment pipeline:
//! 1. Progressive alignment following a guide tree.
//! 2. Iterative refinement with score-based acceptance.
//! 3. Adding sequences to existing alignments.
//!
//! Ports the C `TreeDependentIteration()` from `tditeration.c`,
//! progressive alignment from `disttbfast.c`/`tbfast.c`, and
//! `addonetip()` from `addfunctions.c`.

pub mod progressive;
mod refinement;
mod engine;
mod add;
pub mod external;

pub use progressive::{progressive_align, progressive_align_with_distmtx, MultipleAlignment, StepTrace};
pub use refinement::{iterative_refine, RefinementParams};
pub use engine::{MafftEngine, AlignmentMode};
pub use add::{add_sequences, add_sequences_keeplength};
