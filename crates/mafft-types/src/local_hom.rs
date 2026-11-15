/// A single local homology region between two sequences.
///
/// Replaces the C `LocalHom` linked-list node with a plain struct.
/// Collections of these are stored in a `Vec` rather than a linked list.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HomologyRegion {
    /// Start position in sequence 1.
    pub start1: i32,
    /// End position in sequence 1.
    pub end1: i32,
    /// Start position in sequence 2.
    pub start2: i32,
    /// End position in sequence 2.
    pub end2: i32,
    /// Optimal alignment score for this region.
    pub opt: f64,
    /// Overlap in amino acids.
    pub overlapaa: i32,
    /// Whether this region was extended.
    pub extended: bool,
    /// Importance weight (used in consistency scoring), seen from
    /// sequence 1.
    pub importance: f64,
    /// Importance weight seen from sequence 2 (C's `rimportance`). Equal to
    /// `importance` except between the two passes of `calcimportance`.
    pub rimportance: f64,
    /// 'k' (keep) or 'h' (homolog) classification.
    pub korh: u8,
    /// Remaining count.
    pub nokori: i32,
}

impl Default for HomologyRegion {
    fn default() -> Self {
        Self {
            start1: 0,
            end1: 0,
            start2: 0,
            end2: 0,
            opt: 0.0,
            overlapaa: 0,
            extended: false,
            importance: 0.0,
            rimportance: 0.0,
            korh: b'h',
            nokori: 0,
        }
    }
}

impl HomologyRegion {
    /// The same region seen from the other sequence: sequence-1 and
    /// sequence-2 coordinates exchanged, `importance` and `rimportance`
    /// exchanged.
    #[inline]
    pub fn swapped(self) -> Self {
        Self {
            start1: self.start2,
            end1: self.end2,
            start2: self.start1,
            end2: self.end1,
            importance: self.rimportance,
            rimportance: self.importance,
            ..self
        }
    }
}

/// Table of pairwise local homology information.
///
/// Replaces the C `LocalHom **localhomtable` (linked lists over the
/// `i < j` half of the `njob x njob` matrix, `readlocalhomtable2_half`)
/// with one `Vec` per unordered pair. Each pair's regions are stored once,
/// oriented from the lower index `i` to the higher index `j`:
/// `start1`/`end1` index sequence `i`, `start2`/`end2` sequence `j`,
/// `importance` is the weight seen from `i` and `rimportance` the weight
/// seen from `j` (C's `importance` / `rimportance` on the same node).
/// [`LocalHomologyTable::get`] presents either orientation on read, so
/// callers see `(i, j)` and `(j, i)` as mirror images without a second
/// copy of every region.
#[derive(Debug, Clone, Default)]
pub struct LocalHomologyTable {
    /// Number of sequences.
    pub nseq: usize,
    /// Regions of pair `(i, j)`, `i < j`, at `pair_index(i, j)`.
    entries: Vec<Vec<HomologyRegion>>,
}

/// The regions of one ordered pair `(i, j)`, read in that orientation.
///
/// A view over the stored `min(i, j) < max(i, j)` entry; when `i > j`
/// every region is returned [`HomologyRegion::swapped`].
#[derive(Debug, Clone, Copy)]
pub struct PairRegions<'a> {
    regions: &'a [HomologyRegion],
    swapped: bool,
}

impl<'a> PairRegions<'a> {
    #[inline]
    pub fn len(&self) -> usize {
        self.regions.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.regions.is_empty()
    }

    /// Region `r`, in the pair's orientation.
    #[inline]
    pub fn get(&self, r: usize) -> HomologyRegion {
        let region = self.regions[r];
        if self.swapped { region.swapped() } else { region }
    }

    /// The regions in stored order, in the pair's orientation.
    #[inline]
    pub fn iter(&self) -> PairRegionsIter<'a> {
        PairRegionsIter { inner: self.regions.iter(), swapped: self.swapped }
    }
}

impl<'a> IntoIterator for PairRegions<'a> {
    type Item = HomologyRegion;
    type IntoIter = PairRegionsIter<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Iterator over [`PairRegions`].
#[derive(Debug, Clone)]
pub struct PairRegionsIter<'a> {
    inner: std::slice::Iter<'a, HomologyRegion>,
    swapped: bool,
}

impl Iterator for PairRegionsIter<'_> {
    type Item = HomologyRegion;

    #[inline]
    fn next(&mut self) -> Option<HomologyRegion> {
        let region = *self.inner.next()?;
        Some(if self.swapped { region.swapped() } else { region })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl ExactSizeIterator for PairRegionsIter<'_> {}

impl LocalHomologyTable {
    pub fn new(nseq: usize) -> Self {
        Self {
            nseq,
            entries: vec![Vec::new(); nseq * nseq.saturating_sub(1) / 2],
        }
    }

    /// Position of the unordered pair `{i, j}` (`i != j`) in `entries`:
    /// row-major over the strict upper triangle.
    #[inline]
    fn pair_index(&self, i: usize, j: usize) -> usize {
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        debug_assert!(lo < hi && hi < self.nseq);
        lo * self.nseq - lo * (lo + 1) / 2 + (hi - lo - 1)
    }

    /// Homology regions between sequences `i` and `j`, oriented from `i`
    /// to `j`. Empty for `i == j`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> PairRegions<'_> {
        if i == j {
            return PairRegions { regions: &[], swapped: false };
        }
        PairRegions { regions: &self.entries[self.pair_index(i, j)], swapped: i > j }
    }

    /// The stored regions of pair `(i, j)` in their stored orientation
    /// (from `min(i, j)` to `max(i, j)`), for in-place updates. `i != j`.
    #[inline]
    pub fn stored_mut(&mut self, i: usize, j: usize) -> &mut Vec<HomologyRegion> {
        let idx = self.pair_index(i, j);
        &mut self.entries[idx]
    }

    /// Add a homology region between sequences `i` and `j`, given in the
    /// `(i, j)` orientation. Stored once; `get(j, i)` returns its mirror
    /// image. A region on the diagonal (`i == j`) is dropped.
    pub fn push(&mut self, i: usize, j: usize, region: HomologyRegion) {
        if i == j {
            debug_assert!(false, "local homology region on the diagonal ({i}, {i})");
            return;
        }
        let idx = self.pair_index(i, j);
        self.entries[idx].push(if i < j { region } else { region.swapped() });
    }

    /// Number of stored regions.
    pub fn region_count(&self) -> usize {
        self.entries.iter().map(Vec::len).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mirror_view_swaps_coordinates_and_importance() {
        let mut table = LocalHomologyTable::new(3);
        table.push(2, 0, HomologyRegion {
            start1: 5, end1: 9, start2: 1, end2: 4,
            opt: 3.0, overlapaa: 4, importance: 0.25, rimportance: 0.75,
            ..Default::default()
        });
        // Stored oriented 0 -> 2.
        let fwd = table.get(0, 2);
        assert_eq!(fwd.len(), 1);
        let r = fwd.get(0);
        assert_eq!((r.start1, r.end1, r.start2, r.end2), (1, 4, 5, 9));
        assert_eq!((r.importance, r.rimportance), (0.75, 0.25));
        // Read back as pushed.
        let rev: Vec<HomologyRegion> = table.get(2, 0).iter().collect();
        assert_eq!((rev[0].start1, rev[0].end1, rev[0].start2, rev[0].end2), (5, 9, 1, 4));
        assert_eq!((rev[0].importance, rev[0].rimportance), (0.25, 0.75));
        assert!(table.get(1, 1).is_empty());
        assert!(table.get(0, 1).is_empty());
        assert_eq!(table.region_count(), 1);
    }

    #[test]
    fn pair_index_covers_the_upper_triangle_once() {
        let table = LocalHomologyTable::new(7);
        let mut seen = vec![false; table.entries.len()];
        for i in 0..7 {
            for j in (i + 1)..7 {
                let idx = table.pair_index(i, j);
                assert_eq!(idx, table.pair_index(j, i));
                assert!(!seen[idx], "pair ({i}, {j}) collides");
                seen[idx] = true;
            }
        }
        assert!(seen.iter().all(|&s| s));
    }
}
