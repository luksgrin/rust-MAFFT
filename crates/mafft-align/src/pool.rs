//! Growth of the per-thread scratch buffers that the DP and importance-
//! matrix code reuse across calls.

/// Make `buf` at least `len` elements long, filling new cells with `fill`.
///
/// The buffers this serves (the `width x width` importance matrix and the
/// traceback matrix of the refinement DP) are the largest allocations of
/// the run, and the alignment width moves by a few columns at every
/// accepted refinement step. Two things keep that from costing memory:
///
/// * When the buffer must grow, the old allocation is released before the
///   new one is made instead of going through `realloc`. A `realloc` that
///   cannot extend in place allocates the new block, copies the old
///   contents into it and only then frees the old block, briefly holding
///   both. The copy is also wasted work: every caller rewrites each cell
///   it reads before reading it, so stale contents are never observed.
/// * The new capacity carries a quarter of headroom, so a run of small
///   widenings reuses one allocation instead of freeing and reallocating
///   hundreds of megabytes at each step. Only the `len` cells the caller
///   fills are ever touched, so the headroom costs address space, not
///   resident memory; the point is to free large blocks rarely, because
///   the macOS allocator keeps freed large blocks resident in its cache
///   for a while, and a burst of them was measured to raise the peak RSS
///   of an L-INS-i run an order of magnitude above its live bytes.
#[inline]
pub(crate) fn grow_pooled<T: Copy>(buf: &mut Vec<T>, len: usize, fill: T) {
    if buf.capacity() < len {
        *buf = Vec::new();
        buf.reserve_exact(len + len / 4);
    }
    if buf.len() < len {
        buf.resize(len, fill);
    }
}
