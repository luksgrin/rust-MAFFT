/// A single local homology region between two sequences.
///
/// Replaces the C `LocalHom` linked-list node with a plain struct.
/// Collections of these are stored in a `Vec` rather than a linked list.
#[derive(Debug, Clone, PartialEq)]
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
    /// Importance weight (used in consistency scoring).
    pub importance: f64,
    /// Reverse importance weight.
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

/// Table of pairwise local homology information.
///
/// Replaces the C `LocalHom **localhomtable` (njob x njob linked lists)
/// with a flat map keyed by sequence pair indices.
#[derive(Debug, Clone, Default)]
pub struct LocalHomologyTable {
    /// Number of sequences.
    pub nseq: usize,
    /// Homology regions for each pair (i, j) where i < j.
    /// Indexed as `regions[i * nseq + j]`.
    entries: Vec<Vec<HomologyRegion>>,
}

impl LocalHomologyTable {
    pub fn new(nseq: usize) -> Self {
        Self {
            nseq,
            entries: vec![Vec::new(); nseq * nseq],
        }
    }

    /// Get homology regions between sequences i and j.
    pub fn get(&self, i: usize, j: usize) -> &[HomologyRegion] {
        &self.entries[i * self.nseq + j]
    }

    /// Get mutable homology regions between sequences i and j.
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut Vec<HomologyRegion> {
        &mut self.entries[i * self.nseq + j]
    }

    /// Add a homology region between sequences i and j.
    pub fn push(&mut self, i: usize, j: usize, region: HomologyRegion) {
        self.entries[i * self.nseq + j].push(region);
    }
}

// ---------------------------------------------------------------------------
// FFI conversion: C LocalHom linked list <-> Rust Vec<HomologyRegion>
// ---------------------------------------------------------------------------

#[cfg(feature = "ffi")]
impl HomologyRegion {
    /// Convert a C `LocalHom` linked list into a Vec of Rust regions.
    ///
    /// # Safety
    /// The pointer must point to a valid `LocalHom` linked list (or be null).
    pub unsafe fn from_c_list(head: *const mafft_sys::LocalHom) -> Vec<Self> {
        let mut regions = Vec::new();
        let mut ptr = head;
        while !ptr.is_null() {
            let c = unsafe { &*ptr };
            regions.push(Self {
                start1: c.start1,
                end1: c.end1,
                start2: c.start2,
                end2: c.end2,
                opt: c.opt,
                overlapaa: c.overlapaa,
                extended: c.extended != 0,
                importance: c.importance,
                rimportance: c.rimportance,
                korh: c.korh as u8,
                nokori: c.nokori,
            });
            ptr = c.next;
        }
        regions
    }
}
