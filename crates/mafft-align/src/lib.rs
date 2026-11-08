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
mod multimtx;
mod fft_align;
mod constraints;
mod constrained_align;
mod msalign;

pub use dp::{Alignment, AlignOp, GapModel, matrix_i32_to_f64};
pub use global::global_align;
pub use local::{local_align, release_pairwise_scratch, LocalAlignment};
pub use genaffine::{genaffine_local_align, GenAffineGapModel};
pub use profile::{
    Profile, profile_align, profile_align_imp, profile_align_imp_with_tiebreak,
    profile_align_imp_with_boundary, profile_align_imp_multimtx, BoundaryFreqs,
    profile_align_imp_view, ImpMtx,
    pairwise_align11, pairwise_align11_ex, align_with_anchors,
    reset_dp_pools, reset_cpmx_memo,
};
pub use multimtx::MultiMtx;
pub use fft_align::{fft_profile_align, fft_total_eff, find_fft_anchors, FftAlignParams, Anchor, NucRows};
pub use constraints::{
    build_local_homology_table, build_homology_table,
    build_homology_table_with_unalign, build_imp_diagonal, build_imp_matrix,
    build_imp_matrix_flat_into, take_imp_pool, return_imp_pool,
    recompute_importance, PairAligner, FASTATHRESHOLD_DEFAULT,
    extract_putlocalhom2_regions, build_seed_homology_table,
    merge_homology_tables, parse_hat3_seed, SeedGroup,
};
pub use constrained_align::{
    constrained_profile_align, partial_profile_align, ConstrainedAlignParams,
};
pub use msalign::{msalignmm, msalignmm_with_boundary};

/// Whether any forensic trace is requested through the environment
/// (`RS_*`, `RDBG*`, `MAFFT_DEBUG*`, `MAFFT_DUMP*`, `RUST_MAFFT_TRACE`,
/// `RUST_IMP_DUMP`). Those traces append to files from inside the DP and
/// refinement code in call order, so the parallel paths (which only ever
/// change *where* independent work runs, never its result) fall back to
/// the sequential order while one is set.
pub fn forensic_trace_active() -> bool {
    std::env::vars_os().any(|(k, _)| {
        k.to_str().is_some_and(|k| {
            k.starts_with("RS_")
                || k.starts_with("RDBG")
                || k.starts_with("MAFFT_DEBUG")
                || k.starts_with("MAFFT_DUMP")
                || k == "RUST_MAFFT_TRACE"
                || k == "RUST_IMP_DUMP"
        })
    })
}

thread_local! {
    /// Set by `dp_cancel_scope`: the profile DP on this thread gives up as
    /// soon as `*ptr != expected`.
    static DP_CANCEL: std::cell::Cell<Option<(*const std::sync::atomic::AtomicU64, u64)>> =
        const { std::cell::Cell::new(None) };
}

/// Let the caller abandon a speculative re-alignment early.
///
/// While the returned guard lives, the profile DP on this thread polls
/// `version` every 32 rows and, once it differs from `expected`, stops and
/// returns a throw-away (well-formed, all-gap) alignment. A caller that
/// installs this must discard the result of any computation during which
/// `version` changed; `version` must only ever increase.
pub fn dp_cancel_scope(version: &std::sync::atomic::AtomicU64, expected: u64) -> DpCancelGuard<'_> {
    let prev = DP_CANCEL.with(|c| c.replace(Some((version as *const _, expected))));
    DpCancelGuard { prev, _version: std::marker::PhantomData }
}

/// Restores the previous cancellation token when dropped.
pub struct DpCancelGuard<'a> {
    prev: Option<(*const std::sync::atomic::AtomicU64, u64)>,
    _version: std::marker::PhantomData<&'a std::sync::atomic::AtomicU64>,
}

impl Drop for DpCancelGuard<'_> {
    fn drop(&mut self) {
        DP_CANCEL.with(|c| c.set(self.prev));
    }
}

#[inline]
pub(crate) fn dp_cancelled() -> bool {
    DP_CANCEL.with(|c| match c.get() {
        // SAFETY: only `dp_cancel_scope` installs the pointer, and its guard
        // borrows the atomic for as long as the pointer is installed.
        Some((p, v)) => unsafe { &*p }.load(std::sync::atomic::Ordering::Relaxed) != v,
        None => false,
    })
}

/// Run independent work on the current rayon pool when it has more than
/// one thread and no forensic trace is active.
pub fn parallel_enabled() -> bool {
    rayon::current_num_threads() > 1 && !forensic_trace_active()
}
