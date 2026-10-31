//! Cross-platform deterministic port of FreeBSD/macOS `qsort` (Bentley-
//! McIlroy three-way partitioning quicksort), used to keep tie-break
//! ordering byte-identical to C MAFFT 7.526 *as compiled on macOS*
//! regardless of the host platform the Rust binary runs on.
//!
//! Why this exists: `libc::qsort` delegates to the host C library. macOS
//! uses BSD qsort, Linux uses glibc qsort (merge-sort variant), Windows
//! uses MSVC's qsort. For *truly tied* elements (same score / selfscore /
//! orilen — only when the input contains exactly-duplicate sequences),
//! these implementations disagree on the final order. That makes the
//! Rust binary's output platform-dependent, which fails our parity tests
//! on Linux CI.
//!
//! Algorithm reference: Apple Libc `stdlib/FreeBSD/qsort.c` and
//! `heapsort.c` (in turn from Bentley & McIlroy, "Engineering a Sort
//! Function", 1993). Faithfully preserves swap order so tied elements move
//! identically to the macOS C library.

use std::cmp::Ordering;

/// Sort `slice` in place exactly as macOS libc `qsort` does. The closure
/// `cmp` returns the standard `Ordering` (Less, Equal, Greater) used by
/// Rust's `sort_by`.
///
/// Port of Apple Libc `stdlib/FreeBSD/qsort.c` (`_qsort`, `med3`, `_isort`)
/// including its recursion-depth cap `DEPTH(n) = 2 * (flsl(n) - 1)` that
/// falls back to `heapsort` (`stdlib/FreeBSD/heapsort.c`), the bounded
/// insertion-sort attempt after a swap-free partition, and the
/// smaller-partition-first recursion. Tied elements end up in the same
/// positions as in the C library.
pub fn bsd_qsort<T, F>(slice: &mut [T], mut cmp: F)
where
    T: Clone,
    F: FnMut(&T, &T) -> Ordering,
{
    let n = slice.len();
    let fls = (usize::BITS - n.leading_zeros()) as i64;
    let depth = 2 * (fls - 1);
    qsort_impl(slice, 0, n, &mut cmp, depth);
}

/// `med3`: `a < b ? (b < c ? b : (a < c ? c : a)) : (b > c ? b : (a < c ? a : c))`.
fn med3<T, F>(slice: &[T], a: usize, b: usize, c: usize, cmp: &mut F) -> usize
where
    F: FnMut(&T, &T) -> Ordering,
{
    if cmp(&slice[a], &slice[b]) == Ordering::Less {
        if cmp(&slice[b], &slice[c]) == Ordering::Less {
            b
        } else if cmp(&slice[a], &slice[c]) == Ordering::Less {
            c
        } else {
            a
        }
    } else if cmp(&slice[b], &slice[c]) == Ordering::Greater {
        b
    } else if cmp(&slice[a], &slice[c]) == Ordering::Less {
        a
    } else {
        c
    }
}

/// `vecswap`: swap `n` elements at `a..` with `n` elements at `b..`, front
/// to back.
fn vecswap<T>(slice: &mut [T], a: usize, b: usize, n: usize) {
    for i in 0..n {
        slice.swap(a + i, b + i);
    }
}

/// `_isort`: insertion sort of `slice[start..start + n]`. With a nonzero
/// `swap_limit` it gives up (returning `false`, the range left partially
/// sorted) once more than `swap_limit` swaps were needed.
fn insertion_sort<T, F>(slice: &mut [T], start: usize, n: usize, cmp: &mut F, swap_limit: usize) -> bool
where
    F: FnMut(&T, &T) -> Ordering,
{
    let mut swap_cnt = 0usize;
    for pm in (start + 1)..(start + n) {
        let mut pl = pm;
        while pl > start && cmp(&slice[pl - 1], &slice[pl]) == Ordering::Greater {
            slice.swap(pl, pl - 1);
            swap_cnt += 1;
            if swap_limit != 0 && swap_cnt > swap_limit {
                return false;
            }
            pl -= 1;
        }
    }
    true
}

/// FreeBSD `heapsort` on `s` (1-based heap over `s[0..]`).
fn heapsort<T, F>(s: &mut [T], cmp: &mut F)
where
    T: Clone,
    F: FnMut(&T, &T) -> Ordering,
{
    let n = s.len();
    if n <= 1 {
        return;
    }
    // CREATE
    for l in (1..=n / 2).rev() {
        let mut par_i = l;
        loop {
            let mut child_i = par_i * 2;
            if child_i > n {
                break;
            }
            if child_i < n && cmp(&s[child_i - 1], &s[child_i]) == Ordering::Less {
                child_i += 1;
            }
            if cmp(&s[child_i - 1], &s[par_i - 1]) != Ordering::Greater {
                break;
            }
            s.swap(par_i - 1, child_i - 1);
            par_i = child_i;
        }
    }
    let mut nmemb = n;
    while nmemb > 1 {
        let k = s[nmemb - 1].clone();
        s[nmemb - 1] = s[0].clone();
        nmemb -= 1;
        // SELECT
        let mut par_i = 1usize;
        loop {
            let mut child_i = par_i * 2;
            if child_i > nmemb {
                break;
            }
            if child_i < nmemb && cmp(&s[child_i - 1], &s[child_i]) == Ordering::Less {
                child_i += 1;
            }
            s[par_i - 1] = s[child_i - 1].clone();
            par_i = child_i;
        }
        loop {
            let child_i = par_i;
            par_i = child_i / 2;
            if child_i == 1 || cmp(&k, &s[par_i - 1]) == Ordering::Less {
                s[child_i - 1] = k.clone();
                break;
            }
            s[child_i - 1] = s[par_i - 1].clone();
        }
    }
}

fn qsort_impl<T, F>(slice: &mut [T], mut a: usize, mut n: usize, cmp: &mut F, mut depth_limit: i64)
where
    T: Clone,
    F: FnMut(&T, &T) -> Ordering,
{
    loop {
        let mut swap_cnt = false;
        let limit_hit = depth_limit <= 0;
        depth_limit -= 1;
        if limit_hit {
            heapsort(&mut slice[a..a + n], cmp);
            return;
        }
        if n <= 7 {
            insertion_sort(slice, a, n, cmp, 0);
            return;
        }

        let mut pl = a;
        let mut pm = a + n / 2;
        let mut pn = a + n - 1;
        if n > 40 {
            let d = n / 8;
            pl = med3(slice, pl, pl + d, pl + 2 * d, cmp);
            pm = med3(slice, pm - d, pm, pm + d, cmp);
            pn = med3(slice, pn - 2 * d, pn - d, pn, cmp);
        }
        pm = med3(slice, pl, pm, pn, cmp);
        slice.swap(a, pm);

        // Signed cursors: `pc` may step to `a` (below `pb`).
        let mut pa = a as isize + 1;
        let mut pb = pa;
        let mut pc = (a + n - 1) as isize;
        let mut pd = pc;
        loop {
            while pb <= pc {
                let r = cmp(&slice[pb as usize], &slice[a]);
                if r == Ordering::Greater {
                    break;
                }
                if r == Ordering::Equal {
                    swap_cnt = true;
                    slice.swap(pa as usize, pb as usize);
                    pa += 1;
                }
                pb += 1;
            }
            while pb <= pc {
                let r = cmp(&slice[pc as usize], &slice[a]);
                if r == Ordering::Less {
                    break;
                }
                if r == Ordering::Equal {
                    swap_cnt = true;
                    slice.swap(pc as usize, pd as usize);
                    pd -= 1;
                }
                pc -= 1;
            }
            if pb > pc {
                break;
            }
            slice.swap(pb as usize, pc as usize);
            swap_cnt = true;
            pb += 1;
            pc -= 1;
        }

        let pn_end = (a + n) as isize;
        let ai = a as isize;
        let d1 = (pa - ai).min(pb - pa);
        vecswap(slice, a, (pb - d1) as usize, d1 as usize);
        let d1 = (pd - pc).min(pn_end - pd - 1);
        vecswap(slice, pb as usize, (pn_end - d1) as usize, d1 as usize);

        if !swap_cnt {
            let r = 1 + n / 4;
            if insertion_sort(slice, a, n, cmp, r) {
                return;
            }
        }

        let d1 = (pb - pa) as usize;
        let d2 = (pd - pc) as usize;
        let pn_end = pn_end as usize;
        if d1 <= d2 {
            if d1 > 1 {
                qsort_impl(slice, a, d1, cmp, depth_limit);
            }
            if d2 > 1 {
                a = pn_end - d2;
                n = d2;
                continue;
            }
        } else {
            if d2 > 1 {
                qsort_impl(slice, pn_end - d2, d2, cmp, depth_limit);
            }
            if d1 > 1 {
                n = d1;
                continue;
            }
        }
        return;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_simple_ints() {
        let mut v = vec![5, 3, 8, 1, 4, 7, 2, 6];
        bsd_qsort(&mut v, |a, b| a.cmp(b));
        assert_eq!(v, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn sorts_with_duplicates() {
        let mut v = vec![3, 1, 2, 3, 1, 2, 3, 1];
        bsd_qsort(&mut v, |a, b| a.cmp(b));
        // All sorted, exact tie-break order doesn't matter for primitive ints
        // since equal values are indistinguishable.
        for i in 1..v.len() {
            assert!(v[i - 1] <= v[i]);
        }
    }

    #[test]
    fn handles_small_input() {
        for n in 0..=6 {
            let mut v: Vec<i32> = (0..n).rev().collect();
            bsd_qsort(&mut v, |a, b| a.cmp(b));
            assert_eq!(v, (0..n).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn matches_macos_qsort_on_many_ties() {
        // 60 elements over 5 keys (n > 40: ninther pivot, several
        // partitions of equal elements). Expected tag order recorded from
        // macOS libc `qsort` on the same `(key, tag)` array.
        let mut v: Vec<(usize, usize)> = (0..60).map(|i| ((i * 7) % 5, i)).collect();
        bsd_qsort(&mut v, |a, b| a.0.cmp(&b.0));
        let tags: Vec<usize> = v.iter().map(|e| e.1).collect();
        assert_eq!(tags, vec![
            30, 35, 50, 15, 5, 55, 20, 40, 0, 45, 10, 25, 18, 48, 38, 8, 43, 13, 3, 53,
            23, 28, 33, 58, 46, 56, 26, 51, 41, 6, 11, 1, 36, 21, 31, 16, 59, 4, 9, 14,
            19, 24, 29, 34, 39, 44, 49, 54, 47, 7, 32, 27, 52, 22, 2, 17, 12, 57, 42, 37,
        ]);
    }

    #[test]
    fn matches_macos_qsort_on_partial_tie_pattern() {
        // Mimic the parttree CALL 1 input layout: 36 entries, almost all
        // distinct primary keys, with two entries at indices 0 and 34
        // tied. BSD qsort's median-of-3 + 3-way partition swaps the tied
        // pair so the originally-index-34 element ends up first.
        //
        // Each element is `(key, tag)`; `key` is the primary sort key,
        // `tag` is the "numinseq"-equivalent identity we use to verify
        // post-sort ordering.
        let mut v: Vec<(f64, usize)> = (0..36)
            .map(|i| if i == 0 || i == 34 { (0.0, i) } else { (i as f64, i) })
            .collect();
        // Replicate pick_reference swap: scan for max selfscore, swap to
        // index 0. Here we pre-set the tied entries; just sort.
        bsd_qsort(&mut v, |a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
        // Both tied entries should be at the front.
        assert_eq!(v[0].0, 0.0);
        assert_eq!(v[1].0, 0.0);
        // BSD qsort's specific partition order: the originally-index-34
        // element comes first, originally-index-0 second.
        assert_eq!(v[0].1, 34, "BSD qsort puts originally-index-34 first");
        assert_eq!(v[1].1, 0, "BSD qsort puts originally-index-0 second");
    }
}
