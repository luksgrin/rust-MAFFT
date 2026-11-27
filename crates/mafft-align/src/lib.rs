//! Pairwise and multi-sequence alignment algorithms for MAFFT.
//!
//! Phase 4a: Single-sequence pairwise alignment (Galign11, Lalign11, genalign11).
//! Phase 5: Profile alignment (MSalignmm), FFT-accelerated alignment (Falign),
//!          local homology constraint building (pairlocalalign),
//!          and constrained alignment (Falign_localhom, partSalignmm).

mod dp;
mod pool;
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
pub use profile::dp_critical_path_pays;
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
pub use fft_align::{alignable_site_scores, fft_profile_align, fft_total_eff, find_fft_anchors, FftAlignParams, Anchor, NucRows};
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

thread_local! {
    /// Set by `dp_priority_scope`: `(commit, branch)`. The DP on this
    /// thread spreads its tiles over idle pool threads only while
    /// `*commit == branch`, i.e. while its result is the next one the
    /// refinement walk will decide on.
    static DP_PRIORITY: std::cell::Cell<Option<(*const std::sync::atomic::AtomicUsize, usize)>> =
        const { std::cell::Cell::new(None) };
}

/// Mark the evaluation on this thread as speculative branch `branch` of a
/// refinement cycle whose commit point is `*commit`.
///
/// While the guard lives, [`dp_inner_parallel`] is true only when
/// `*commit == branch`: the branch at the commit point is on the critical
/// path of the cycle and may recruit idle pool threads for its DP; a
/// branch further ahead is speculative and runs on its own thread, so the
/// idle threads stay available to the critical one. `commit` must only
/// ever increase. The decision affects scheduling only, never a result.
pub fn dp_priority_scope(
    commit: &std::sync::atomic::AtomicUsize,
    branch: usize,
) -> DpPriorityGuard<'_> {
    let prev = DP_PRIORITY.with(|c| c.replace(Some((commit as *const _, branch))));
    DpPriorityGuard { prev, _commit: std::marker::PhantomData }
}

/// Restores the previous priority token when dropped.
pub struct DpPriorityGuard<'a> {
    prev: Option<(*const std::sync::atomic::AtomicUsize, usize)>,
    _commit: std::marker::PhantomData<&'a std::sync::atomic::AtomicUsize>,
}

impl Drop for DpPriorityGuard<'_> {
    fn drop(&mut self) {
        DP_PRIORITY.with(|c| c.set(self.prev));
    }
}

/// May the DP running on this thread hand work to idle pool threads?
///
/// True on a multi-thread pool when no [`dp_priority_scope`] is installed
/// (the progressive stage, segment refinement) or when the installed one
/// says this evaluation is at the commit point. Cheap enough to poll at
/// block boundaries, so an evaluation that becomes critical part-way
/// through picks up the idle threads from its next block.
#[inline]
pub fn dp_inner_parallel() -> bool {
    parallel_enabled()
        && DP_PRIORITY.with(|c| match c.get() {
            // SAFETY: only `dp_priority_scope` installs the pointer, and its
            // guard borrows the atomic for as long as the pointer is installed.
            Some((p, k)) => unsafe { &*p }.load(std::sync::atomic::Ordering::Acquire) == k,
            None => true,
        })
}

/// This thread's cancellation token, for work this thread hands to other
/// pool threads (the tiles of `profile::dp_fill_tiled`), which have no
/// token of their own. Valid while this thread's `DpCancelGuard` lives,
/// i.e. for the duration of the DP call that spawned the work.
#[derive(Clone, Copy)]
pub(crate) struct DpCancelToken(Option<(*const std::sync::atomic::AtomicU64, u64)>);
// SAFETY: the pointee is an `AtomicU64` that outlives every holder of the
// token (see `dp_cancel_token`); atomics are safe to read from any thread.
unsafe impl Send for DpCancelToken {}
unsafe impl Sync for DpCancelToken {}

impl DpCancelToken {
    #[inline]
    pub(crate) fn cancelled(self) -> bool {
        match self.0 {
            // SAFETY: see `dp_cancelled`.
            Some((p, v)) => unsafe { &*p }.load(std::sync::atomic::Ordering::Relaxed) != v,
            None => false,
        }
    }
}

#[inline]
pub(crate) fn dp_cancel_token() -> DpCancelToken {
    DpCancelToken(DP_CANCEL.with(|c| c.get()))
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
