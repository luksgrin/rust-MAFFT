//! Growth of the per-thread scratch buffers that the DP and importance-
//! matrix code reuse across calls.

/// Make `buf` at least `len` elements long, filling new cells with `fill`.
///
/// The buffers this serves (the `width x width` importance matrix and the
/// traceback matrices of the refinement and progressive DPs) are the
/// largest allocations of the run. Two access patterns meet here: in
/// refinement the alignment width moves by a few columns at every accepted
/// step, in the progressive stage the DP grows monotonically towards the
/// final width. The capacity therefore carries a quarter of headroom: a
/// run of small widenings reuses one allocation instead of touching a
/// block of hundreds of megabytes at each step, and a monotone growth
/// re-allocates O(log) times instead of once per merge. Only the `len`
/// cells the caller fills are ever touched, so the headroom costs address
/// space, not resident memory.
///
/// Growth goes through `reserve_exact`, i.e. `realloc`: on macOS a large
/// block is usually extended in place, and when it is not the old block is
/// handed back to the allocator with the new one already in use, which was
/// measured to leave less resident than freeing first and allocating a
/// fresh block (the freed block is kept resident in the allocator's
/// large-block cache, so a sequence of grow steps that frees first parks
/// one dead block per step; PF00005 n100 `--parttree` peaked at 116 MB
/// that way against 48 MB with `realloc`). Every caller rewrites each cell
/// it reads before reading it, so the preserved contents are never
/// observed either way.
#[inline]
pub(crate) fn grow_pooled<T: Copy>(buf: &mut Vec<T>, len: usize, fill: T) {
    if buf.capacity() < len {
        buf.reserve_exact(len + len / 4 - buf.len());
    }
    if buf.len() < len {
        buf.resize(len, fill);
    }
}
