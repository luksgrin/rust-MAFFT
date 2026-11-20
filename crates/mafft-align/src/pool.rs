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

/// Distance between the page offsets [`coloured_window`] hands out.
const COLOUR_STEP: usize = 512;
const PAGE: usize = 4096;

/// Grow `buf` to hold a window of `len` cells whose first cell lies
/// `colour * 512` bytes past a 4 KiB boundary, and return the window's
/// start index (`buf[start..start + len]`). Cells of a grown buffer are
/// filled with `fill`; cells already there keep their contents, so the
/// caller writes every window cell it reads.
///
/// The profile DP's inner loop (`profile.rs`, `j_loop`) walks eight row
/// buffers in lock step: column `j` stores to `currentw`, `mj`, `mpj` and
/// the row of `h`, and the next column loads `previousw[j]`, `fgcp2[j]`
/// and `nongap_freq[j]`. The rows are heap blocks of up to ~90 KB that
/// the allocator returns page-aligned and, within one size class, at
/// regular spacings, so which of them shared their low address bits was
/// decided by the history of the process heap. With `mj` and `fgcp2`
/// exactly 2 MiB apart, the load of `fgcp2[j]` right after the store to
/// `mj[j]` matched it in the low 21 address bits, and the loop took 42 %
/// more cycles for the same instructions (PF00069 n=250, FFT-NS-2: 8.8 s
/// against 5.9 s); with every row on the same 16 KiB offset they also
/// compete for the same L1 sets. Giving each row its own offset in the
/// page (colours 512 bytes apart) keeps the low 12 bits of any store and
/// load of the loop at least 504 bytes apart, whatever the allocator
/// returns. It costs at most 4 KiB of slack per buffer.
#[inline]
pub(crate) fn coloured_window<T: Copy>(
    buf: &mut Vec<T>,
    len: usize,
    fill: T,
    colour: usize,
) -> usize {
    let size = std::mem::size_of::<T>();
    debug_assert!(size.is_power_of_two() && size <= COLOUR_STEP);
    grow_pooled(buf, len + PAGE / size, fill);
    let want = (colour * COLOUR_STEP) % PAGE;
    let have = buf.as_ptr() as usize % PAGE;
    // Both offsets are multiples of `size` (the allocation is aligned to
    // `T`), so the distance is a whole number of cells, below PAGE / size.
    ((want + PAGE - have) % PAGE) / size
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coloured_window_lands_on_its_page_offset() {
        let mut f: Vec<f64> = Vec::new();
        let mut u: Vec<usize> = Vec::new();
        for (len, colour) in [(1usize, 0usize), (953, 3), (8885, 7), (100, 9), (20000, 5)] {
            let at = coloured_window(&mut f, len, 0.0, colour);
            assert!(at + len <= f.len());
            assert_eq!(f[at..].as_ptr() as usize % PAGE, colour * COLOUR_STEP % PAGE);
            let at = coloured_window(&mut u, len, 0, colour);
            assert!(at + len <= u.len());
            assert_eq!(u[at..].as_ptr() as usize % PAGE, colour * COLOUR_STEP % PAGE);
        }
    }
}
