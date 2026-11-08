/// Pairwise distance computation between sequences.

/// A symmetric distance matrix stored as upper triangle.
///
/// `data[i]` has `nseq - i - 1` elements, storing distances from
/// sequence `i` to sequences `i+1, i+2, ..., nseq-1`.
#[derive(Debug, Clone)]
pub struct DistanceMatrix {
    pub nseq: usize,
    data: Vec<Vec<f64>>,
}

impl DistanceMatrix {
    pub fn new(nseq: usize) -> Self {
        let data = (0..nseq)
            .map(|i| vec![0.0; nseq - i - 1])
            .collect();
        Self { nseq, data }
    }

    /// Create from a full symmetric matrix.
    pub fn from_full(matrix: &[Vec<f64>]) -> Self {
        let nseq = matrix.len();
        let mut dm = Self::new(nseq);
        for i in 0..nseq {
            for j in (i + 1)..nseq {
                dm.set(i, j, matrix[i][j]);
            }
        }
        dm
    }

    /// Get distance between sequences i and j.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            0.0
        } else if i < j {
            self.data[i][j - i - 1]
        } else {
            self.data[j][i - j - 1]
        }
    }

    /// Set distance between sequences i and j (symmetric).
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        if i < j {
            self.data[i][j - i - 1] = val;
        } else if j < i {
            self.data[j][i - j - 1] = val;
        }
    }

    /// Build from the upper triangle given row by row: `rows[i][j - i - 1]`
    /// is the distance between `i` and `j > i` (the internal layout).
    pub fn from_upper_rows(rows: Vec<Vec<f64>>) -> Self {
        let nseq = rows.len();
        for (i, r) in rows.iter().enumerate() {
            assert_eq!(r.len(), nseq - i - 1, "row {i} of the upper triangle has the wrong length");
        }
        Self { nseq, data: rows }
    }

}

/// Compute identity-based distance between two aligned sequences.
///
/// Distance = 1.0 - (identical positions / aligned positions excluding gaps).
///
/// Uses branchless counting for auto-vectorization: the gap check and
/// equality check are converted to integer masks that LLVM can vectorize
/// with SIMD compare+accumulate instructions.
pub fn pairwise_identity_distance(seq1: &[u8], seq2: &[u8]) -> f64 {
    let mut matches = 0u64;
    let mut aligned = 0u64;

    let len = seq1.len().min(seq2.len());
    for k in 0..len {
        let a = seq1[k];
        let b = seq2[k];
        // Branchless: both_nongap = 1 if neither is '-', else 0
        let both_nongap = ((a != b'-') & (b != b'-')) as u64;
        // Branchless: is_match = 1 if a == b, else 0
        let is_match = (a == b) as u64;
        aligned += both_nongap;
        matches += both_nongap & is_match;
    }

    if aligned == 0 {
        1.0
    } else {
        1.0 - matches as f64 / aligned as f64
    }
}

/// Compute the k-tuple (6-mer) distance between two unaligned sequences
/// whose alphabet is already known.
///
/// Ports the C `commonsextet_p` + `distcompact` algorithm exactly:
/// - Amino acids are mapped to 6 physicochemical groups (matching C's `amino_grp`)
/// - DNA bases are mapped to 4 groups; `n` and the IUPAC ambiguity codes
///   are dropped from the 6-mer stream (`disttbfast.c::seq_grp_nuc`,
///   `if( tmp < 4 )`, with `locgrpn` from `DNA.h:44`)
/// - K-tuples are encoded as integers into a flat frequency table
/// - Distance includes C's length adjustment factor (`lenfac`)
///
/// Formula: `(1.0 - common/min(ss1,ss2)) * lenfac * 2.0`
///
/// `is_nucleotide` is the input-wide sequence type (`dorp` in C, decided
/// once per run from the whole input by `countlen`/`getnumlen`), never a
/// per-pair guess: C never looks at the residues of a pair to choose the
/// alphabet. Production callers pass `ScoringContext::seq_type.is_nucleotide()`.
pub fn ktuple_distance_typed(seq1: &[u8], seq2: &[u8], k: usize, is_nucleotide: bool) -> f64 {
    if seq1.is_empty() || seq2.is_empty() || k == 0 {
        return 2.0;
    }
    if seq1.len() < k || seq2.len() < k {
        return 2.0;
    }
    if is_nucleotide {
        ktuple_distance_nuc(seq1, seq2, k)
    } else {
        ktuple_distance_aa(seq1, seq2, k)
    }
}

/// [`ktuple_distance_typed`] with the alphabet sniffed from the first 100
/// letters of the pair.
///
/// Only for callers that have no input-wide sequence type at hand (tests,
/// diagnostics). The sniff is per pair, so a nucleotide sequence with an
/// IUPAC ambiguity code (R, Y, K, M, W, S, B, D, H, V) among its first 100
/// bases is treated as protein, and every pair it leads gets the protein
/// 6-group distance: with the protein grouping A, G and T all fall into
/// group 0, the pair looks far closer than it is, and the guide tree moves.
/// C never does this (see [`ktuple_distance_typed`]). The engine and the
/// `--add` path therefore call the typed function with the run-wide
/// `SeqType`; this wrapper must not be used on production paths.
pub fn ktuple_distance(seq1: &[u8], seq2: &[u8], k: usize) -> f64 {
    let is_dna = seq1.iter().chain(seq2.iter())
        .filter(|&&c| c.is_ascii_alphabetic())
        .take(100)
        .all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u' | b'N' | b'n'));
    ktuple_distance_typed(seq1, seq2, k, is_dna)
}

/// Amino acid group mapping (C's `locgrpd` from blosum.c).
///
/// Maps each amino acid to one of 6 physicochemical groups:
///   0: A, G, P, S, T (small/turn)
///   1: I, L, M, V, J (hydrophobic)
///   2: N, D, Q, E, B, Z (charged/polar)
///   3: R, H, K (positive)
///   4: F, W, Y (aromatic)
///   5: C (cysteine)
/// Values >= 6 are skipped (X, gaps, etc.)
fn amino_group(c: u8) -> Option<u32> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'R' => Some(3),
        b'N' => Some(2),
        b'D' => Some(2),
        b'C' => Some(5),
        b'Q' => Some(2),
        b'E' => Some(2),
        b'G' => Some(0),
        b'H' => Some(3),
        b'I' => Some(1),
        b'L' => Some(1),
        b'K' => Some(3),
        b'M' => Some(1),
        b'F' => Some(4),
        b'P' => Some(0),
        b'S' => Some(0),
        b'T' => Some(0),
        b'W' => Some(4),
        b'Y' => Some(4),
        b'V' => Some(1),
        b'B' => Some(2),
        b'Z' => Some(2),
        b'J' => Some(1),
        _ => None, // X, gaps, unknown → skip
    }
}

/// Scoring-matrix-based pairwise distance from aligned sequences.
///
/// Ports C's `naivepairscore11` + distance formula. Used for guide tree
/// rebuilding in the retree pass (instead of simple identity distance).
///
/// Algorithm: strip common gap columns, then score aligned positions using
/// the substitution matrix. Gap blocks get a flat penalty (once per block).
///
/// Distance = (1.0 - score / min(selfscore_i, selfscore_j)) * 2.0
pub fn scoring_matrix_distance(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    penalty: i32,
) -> f64 {
    let ss1 = self_score(seq1, matrix, amino_map, penalty);
    let ss2 = self_score(seq2, matrix, amino_map, penalty);
    let bunbo = if ss1 < ss2 { ss1 } else { ss2 };
    if bunbo == 0.0 {
        return 2.0;
    }
    let score = naive_pair_score(seq1, seq2, matrix, amino_map, penalty);
    // C's distcompact_msa clamps high at 10.0 and low at 0.0 (with warning).
    let mut d = (1.0 - score / bunbo) * 2.0;
    if d > 10.0 { d = 10.0; }
    if d < 0.0 { d = 0.0; }
    d
}


/// C's naivepairscore11: score two aligned sequences.
fn naive_pair_score(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    penalty: i32,
) -> f64 {
    let nalpha = matrix.len();
    // Strip common gap columns
    let mut s1 = Vec::new();
    let mut s2 = Vec::new();
    let len = seq1.len().min(seq2.len());
    for k in 0..len {
        if seq1[k] != b'-' || seq2[k] != b'-' {
            s1.push(seq1[k]);
            s2.push(seq2[k]);
        }
    }

    let mut score = 0.0f64;
    let mut k = 0;
    while k < s1.len() {
        if s1[k] == b'-' {
            score += penalty as f64;
            while k < s1.len() && s1[k] == b'-' { k += 1; }
            continue;
        }
        if s2[k] == b'-' {
            score += penalty as f64;
            while k < s2.len() && s2[k] == b'-' { k += 1; }
            continue;
        }
        let i = amino_map[s1[k] as usize] as usize;
        let j = amino_map[s2[k] as usize] as usize;
        if i < nalpha && j < nalpha {
            score += matrix[i][j] as f64;
        }
        k += 1;
    }
    score
}

/// One aligned row prepared for [`MsaDistance::distance`]: its self-score
/// and the columns that hold a residue, computed once per row.
///
/// [`scoring_matrix_distance`] recomputes both self-scores and copies the
/// two rows without their common-gap columns for every pair, so a matrix
/// costs about three passes over the full alignment width per pair. C
/// computes `selfscore[i]` once per sequence (`disttbfast.c:4454`,
/// `dndpre.c:412`) and walks only the non-common-gap columns through
/// `makeskiptable` + `commongappickpairfast` (`mltaln9.c:15243`,
/// `mltaln9.c:13322`).
#[derive(Debug, Clone)]
pub struct ScoredRow {
    len: usize,
    self_score: f64,
    /// Columns that are not `'-'`, ascending.
    pos: Vec<u32>,
    /// `amino_map` of the residue in each of those columns.
    code: Vec<u8>,
}

/// [`scoring_matrix_distance`] for many pairs of rows of one alignment.
///
/// Returns bit-identical values: the pair score adds the same terms in the
/// same column order (a match score per column where both rows have a
/// residue, `penalty` once per run of columns where exactly one row has a
/// gap, columns where both have a gap skipped), and the self-scores come
/// from the same function, called once per row instead of once per pair.
pub struct MsaDistance<'m> {
    matrix: &'m [Vec<i32>],
    amino_map: &'m [u8; 256],
    penalty: i32,
    nalpha: usize,
    /// `matrix` flattened row-major (`nalpha * nalpha`).
    flat: Vec<i32>,
}

impl<'m> MsaDistance<'m> {
    /// `None` when `matrix` is not square (the per-pair function is then
    /// the only safe choice, since it indexes `matrix[i][j]` lazily).
    pub fn new(matrix: &'m [Vec<i32>], amino_map: &'m [u8; 256], penalty: i32) -> Option<Self> {
        let nalpha = matrix.len();
        if matrix.iter().any(|row| row.len() != nalpha) {
            return None;
        }
        let flat = matrix.iter().flat_map(|row| row.iter().copied()).collect();
        Some(Self { matrix, amino_map, penalty, nalpha, flat })
    }

    /// Prepare one aligned row.
    pub fn row(&self, seq: &[u8]) -> ScoredRow {
        assert!(seq.len() <= u32::MAX as usize, "alignment too wide");
        let mut pos = Vec::with_capacity(seq.len() / 2);
        let mut code = Vec::with_capacity(seq.len() / 2);
        for (k, &c) in seq.iter().enumerate() {
            if c != b'-' {
                pos.push(k as u32);
                code.push(self.amino_map[c as usize]);
            }
        }
        ScoredRow {
            len: seq.len(),
            self_score: self_score(seq, self.matrix, self.amino_map, self.penalty),
            pos,
            code,
        }
    }

    /// `scoring_matrix_distance(seq_a, seq_b, ..)` for prepared rows.
    pub fn distance(&self, a: &ScoredRow, b: &ScoredRow) -> f64 {
        let ss1 = a.self_score;
        let ss2 = b.self_score;
        let bunbo = if ss1 < ss2 { ss1 } else { ss2 };
        if bunbo == 0.0 {
            return 2.0;
        }
        let score = self.pair_score(a, b);
        let mut d = (1.0 - score / bunbo) * 2.0;
        if d > 10.0 { d = 10.0; }
        if d < 0.0 { d = 0.0; }
        d
    }

    /// `naive_pair_score` without materialising the common-gap-stripped
    /// rows: merge the two residue-column lists. A column in both lists is
    /// a match; a column in only one list is a gap in the other row, and
    /// consecutive such columns of the same kind form one penalised run.
    fn pair_score(&self, a: &ScoredRow, b: &ScoredRow) -> f64 {
        const MATCH: u8 = 0;
        const GAP1: u8 = 1; // row a has the gap
        const GAP2: u8 = 2; // row b has the gap
        // `naive_pair_score` only looks at the first min(len) columns.
        let lim = a.len.min(b.len) as u32;
        let na = a.pos.partition_point(|&p| p < lim);
        let nb = b.pos.partition_point(|&p| p < lim);
        let (pos1, code1) = (&a.pos[..na], &a.code[..na]);
        let (pos2, code2) = (&b.pos[..nb], &b.code[..nb]);
        let nalpha = self.nalpha;
        let pen = self.penalty as f64;
        let mut score = 0.0f64;
        let mut state = MATCH;
        let (mut i, mut j) = (0usize, 0usize);
        while i < na && j < nb {
            let p1 = pos1[i];
            let p2 = pos2[j];
            if p1 == p2 {
                let c1 = code1[i] as usize;
                let c2 = code2[j] as usize;
                if c1 < nalpha && c2 < nalpha {
                    score += self.flat[c1 * nalpha + c2] as f64;
                }
                state = MATCH;
                i += 1;
                j += 1;
            } else if p1 < p2 {
                if state != GAP2 {
                    score += pen;
                    state = GAP2;
                }
                i += 1;
            } else {
                if state != GAP1 {
                    score += pen;
                    state = GAP1;
                }
                j += 1;
            }
        }
        // Whatever is left of one list is a single run of gaps in the other row.
        if i < na && state != GAP2 {
            score += pen;
        }
        if j < nb && state != GAP1 {
            score += pen;
        }
        score
    }
}

/// One sequence prepared for [`ktuple_distances_from`]: what
/// `ktuple_distance_typed(seq, other, k, is_nucleotide)` derives from `seq`
/// alone.
///
/// `ktuple_distance_typed` rebuilds both sequences' group strings and
/// k-tuple point vectors, and allocates two zeroed `alphabet^k`-entry
/// tables (two 186 KB `calloc`s for protein), for every pair. C builds
/// `pointt[i]` once per sequence (`disttbfast.c:3601-3633`), one
/// composition table per row (`disttbfast.c:3827-3835`), and
/// `commonsextet_p` resets only the `memo` entries it touched
/// (`mltaln9.c:15376`). [`ktuple_distances_from`] does the same.
///
/// The alphabet is the run-wide sequence type (`dorp` in C), fixed once
/// for the whole matrix by the caller; it is never inferred from the
/// residues of a sequence or a pair.
pub struct KtupleSeq {
    k: usize,
    /// `true` for the 4-letter nucleotide alphabet, `false` for the
    /// 6-group protein alphabet.
    dna: bool,
    /// Residue count without `'-'` (`nogaplen` in C), for `lenfac`.
    nogap: usize,
    /// Whether the raw sequence is empty or shorter than `k` bytes (the
    /// early `2.0` returns of `ktuple_distance_typed`).
    too_short: bool,
    /// k-tuple codes (`make_point_table`) in the chosen alphabet, or
    /// `None` when the sequence has fewer than `k` grouped letters.
    points: Option<Vec<u32>>,
}

impl KtupleSeq {
    /// Prepare `seq` for the alphabet given by `is_nucleotide`
    /// (`ScoringContext::seq_type.is_nucleotide()` on production paths).
    pub fn new(seq: &[u8], k: usize, is_nucleotide: bool) -> Self {
        let too_short = seq.is_empty() || k == 0 || seq.len() < k;
        let points = if too_short {
            None
        } else {
            let (group_fn, alphabet): (fn(u8) -> Option<u32>, u32) =
                if is_nucleotide { (nuc_group, 4) } else { (amino_group, 6) };
            let groups = seq_to_groups(seq, group_fn);
            if groups.len() < k {
                None
            } else {
                Some(make_point_table(&groups, k, alphabet))
            }
        };
        Self { k, dna: is_nucleotide, nogap: nogap_len(seq), too_short, points }
    }
}

/// Per-thread scratch for [`ktuple_distances_from`]: the composition table
/// of the current row's sequence (`table1` in C) and the `memo` counts of
/// `commonsextet_p`, both kept all-zero between uses by resetting only the
/// entries that were touched.
#[derive(Default)]
pub struct KtupleScratch {
    table: Vec<u32>,
    /// The row index and the codes counted into `table`.
    loaded: Option<(usize, Vec<u32>)>,
    memo: Vec<u32>,
}

impl KtupleScratch {
    /// `common_sextets(make_composition_table(p1), p2)`: the same counting
    /// loop, with `table1` built once per row instead of once per pair.
    fn common(&mut self, row: usize, dna: bool, k: usize, p1: &[u32], p2: &[u32]) -> u64 {
        let tsize = if dna { 4usize } else { 6usize }.pow(k as u32);
        if self.table.len() < tsize {
            self.table = vec![0; tsize];
            self.loaded = None;
        }
        if self.memo.len() < tsize {
            self.memo.resize(tsize, 0);
        }
        if self.loaded.as_ref().map(|(r, _)| *r) != Some(row) {
            let table = &mut self.table;
            if let Some((_, prev)) = self.loaded.take() {
                for &p in &prev {
                    table[p as usize] = 0;
                }
            }
            for &p in p1 {
                table[p as usize] += 1;
            }
            self.loaded = Some((row, p1.to_vec()));
        }
        let table = &self.table;
        let memo = &mut self.memo;
        let mut value = 0u64;
        for &p in p2 {
            let tmp = memo[p as usize];
            memo[p as usize] = tmp + 1;
            if tmp < table[p as usize] {
                value += 1;
            }
        }
        for &p in p2 {
            memo[p as usize] = 0;
        }
        value
    }
}

/// `ktuple_distance_typed(seqs[i], seqs[j], k, is_nucleotide)` for every
/// `j > i`, bit-identical, from sequences prepared with the same `k` and
/// alphabet.
pub fn ktuple_distances_from(prepared: &[KtupleSeq], i: usize, scratch: &mut KtupleScratch) -> Vec<f64> {
    let a = &prepared[i];
    let k = a.k;
    let dna = a.dna;
    prepared[i + 1..]
        .iter()
        .map(|b| {
            debug_assert_eq!(k, b.k);
            debug_assert_eq!(dna, b.dna);
            if a.too_short || b.too_short {
                return 2.0;
            }
            let (Some(p1), Some(p2)) = (a.points.as_deref(), b.points.as_deref()) else {
                return 2.0;
            };
            let ss1 = p1.len() as u64;
            let ss2 = p2.len() as u64;
            if ss1 == 0 || ss2 == 0 {
                return 2.0;
            }
            let common = scratch.common(i, dna, k, p1, p2);
            let (len1, len2) = (a.nogap as f64, b.nogap as f64);
            let lenfac = if dna {
                compute_lenfac(len1, len2, 0.01, 2500.0, 2500.0, 0.1)
            } else {
                compute_lenfac(len1, len2, 0.01, 10000.0, 10000.0, 0.1)
            };
            (1.0 - common as f64 / ss1.min(ss2) as f64) * lenfac * 2.0
        })
        .collect()
}

/// Self-score: naivepairscore11(seq, seq) = sum of diagonal matrix values.
fn self_score(
    seq: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    _penalty: i32,
) -> f64 {
    let nalpha = matrix.len();
    let mut score = 0.0f64;
    for &ch in seq {
        if ch != b'-' {
            let i = amino_map[ch as usize] as usize;
            if i < nalpha {
                score += matrix[i][i] as f64;
            }
        }
    }
    score
}

/// DNA base group mapping: A→0, C→1, G→2, T/U→3.
fn nuc_group(c: u8) -> Option<u32> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// Convert a sequence to group-encoded integers, skipping unknowns.
fn seq_to_groups(seq: &[u8], group_fn: fn(u8) -> Option<u32>) -> Vec<u32> {
    seq.iter().filter_map(|&c| group_fn(c)).collect()
}

/// Build the point vector (encoded k-tuple indices) from a group sequence.
///
/// For alphabet_size=6, k=6: point = g[0]*7776 + g[1]*1296 + g[2]*216 + g[3]*36 + g[4]*6 + g[5]
fn make_point_table(groups: &[u32], k: usize, alphabet_size: u32) -> Vec<u32> {
    if groups.len() < k {
        return Vec::new();
    }

    let top_power = alphabet_size.pow((k - 1) as u32);
    let mut points = Vec::with_capacity(groups.len() - k + 1);

    // Compute first k-tuple
    let mut point: u32 = 0;
    for i in 0..k {
        point = point * alphabet_size + groups[i];
    }
    points.push(point);

    // Sliding window
    for i in k..groups.len() {
        point = (point - groups[i - k] * top_power) * alphabet_size + groups[i];
        points.push(point);
    }

    points
}

/// Build composition table (frequency of each k-tuple index).
fn make_composition_table(points: &[u32], tsize: usize) -> Vec<u32> {
    let mut table = vec![0u32; tsize];
    for &p in points {
        table[p as usize] += 1;
    }
    table
}

/// Count common sextets: C's `commonsextet_p` algorithm.
///
/// For each k-tuple in seq2's point vector, checks if seq1's composition
/// table has remaining count. Uses min-capping to avoid double-counting.
fn common_sextets(table1: &[u32], points2: &[u32], tsize: usize) -> u64 {
    let mut memo = vec![0u32; tsize];
    let mut value = 0u64;

    for &point in points2 {
        let p = point as usize;
        let tmp = memo[p];
        memo[p] += 1;
        if tmp < table1[p] {
            value += 1;
        }
    }

    value
}

/// Protein k-tuple distance with C's exact algorithm.
fn ktuple_distance_aa(seq1: &[u8], seq2: &[u8], k: usize) -> f64 {
    let groups1 = seq_to_groups(seq1, amino_group);
    let groups2 = seq_to_groups(seq2, amino_group);

    if groups1.len() < k || groups2.len() < k {
        return 2.0;
    }

    let alphabet_size: u32 = 6;
    let tsize = alphabet_size.pow(k as u32) as usize; // 6^6 = 46656

    let points1 = make_point_table(&groups1, k, alphabet_size);
    let points2 = make_point_table(&groups2, k, alphabet_size);

    let table1 = make_composition_table(&points1, tsize);

    let ss1 = points1.len() as u64;
    let ss2 = points2.len() as u64;

    if ss1 == 0 || ss2 == 0 {
        return 2.0;
    }

    let common = common_sextets(&table1, &points2, tsize);

    // C's disttbfast.c:3845-3856 uses `nogaplen` (gaps-only stripped),
    // NOT the filtered group length, for lenfac. Sequences with X / '.' /
    // other unknowns contribute to nogaplen but are skipped from groups.
    let len1 = nogap_len(seq1) as f64;
    let len2 = nogap_len(seq2) as f64;
    let lenfac = compute_lenfac(len1, len2, 0.01, 10000.0, 10000.0, 0.1);

    // C `disttbfast.c:3894` stores `(1 - common/min(ss1,ss2)) * lenfac * 2`
    // without clamping. Distances > 2.0 occur for very-different-length
    // pairs and affect UPGMA tie-breaks. Earlier clamping caused BB40041
    // (mixed-length 105..1084 family) to build a different tree.
    (1.0 - common as f64 / ss1.min(ss2) as f64) * lenfac * 2.0
}

/// DNA k-tuple distance with C's exact algorithm.
fn ktuple_distance_nuc(seq1: &[u8], seq2: &[u8], k: usize) -> f64 {
    let groups1 = seq_to_groups(seq1, nuc_group);
    let groups2 = seq_to_groups(seq2, nuc_group);

    if groups1.len() < k || groups2.len() < k {
        return 2.0;
    }

    let alphabet_size: u32 = 4;
    let tsize = alphabet_size.pow(k as u32) as usize; // 4^6 = 4096

    let points1 = make_point_table(&groups1, k, alphabet_size);
    let points2 = make_point_table(&groups2, k, alphabet_size);

    let table1 = make_composition_table(&points1, tsize);

    let ss1 = points1.len() as u64;
    let ss2 = points2.len() as u64;

    if ss1 == 0 || ss2 == 0 {
        return 2.0;
    }

    let common = common_sextets(&table1, &points2, tsize);

    // C's disttbfast.c:3845-3856 uses `nogaplen` for lenfac (see aa version).
    let len1 = nogap_len(seq1) as f64;
    let len2 = nogap_len(seq2) as f64;
    let lenfac = compute_lenfac(len1, len2, 0.01, 2500.0, 2500.0, 0.1);

    // No clamp — C doesn't cap (see ktuple_distance_aa comment).
    (1.0 - common as f64 / ss1.min(ss2) as f64) * lenfac * 2.0
}

fn nogap_len(seq: &[u8]) -> usize {
    seq.iter().filter(|&&c| c != b'-').count()
}

/// C's length adjustment factor.
///
/// `lenfac = 1.0 / (shorter/longer * d + b/(longer + c) + a)`
/// (`disttbfast.c:3856`). Apple clang at `-O3` contracts
/// `shorter/longer * d + b/(longer + c)` into one `fmadd` and then adds `a`
/// (checked in the generated assembly). The unfused sum is 1 ULP low for
/// some unequal lengths (29 vs 30 with the protein constants); that ULP
/// moves UPGMA branch lengths and flips a progressive tie on large sets of
/// short sequences.
fn compute_lenfac(len1: f64, len2: f64, a: f64, b: f64, c: f64, d: f64) -> f64 {
    let longer = len1.max(len2);
    let shorter = len1.min(len2);
    if longer == 0.0 {
        return 1.0;
    }
    let ratio = shorter / longer;
    let mid = b / (longer + c);
    1.0 / (mafft_types::fp::fmadd(ratio, d, mid) + a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_distance_identical() {
        assert!((pairwise_identity_distance(b"ACGT", b"ACGT")).abs() < 1e-10);
    }

    #[test]
    fn identity_distance_half() {
        // 2 out of 4 match
        assert!((pairwise_identity_distance(b"ACGT", b"ACGA") - 0.25).abs() < 1e-10);
    }

    #[test]
    fn identity_distance_with_gaps() {
        // Gaps excluded: A-GT vs A-GT = 3 aligned, 3 match
        assert!((pairwise_identity_distance(b"A-GT", b"A-GT")).abs() < 1e-10);
    }

    #[test]
    fn ktuple_identical() {
        let d = ktuple_distance(b"ACGTACGTACGT", b"ACGTACGTACGT", 6);
        assert!(d < 0.01, "distance should be ~0 for identical, got {d}");
    }

    #[test]
    fn ktuple_different() {
        // ACDEFGHIKLMN spans many k-tuple groups; WWWWWWWWWWWW is all-aromatic
        // (group 4). No common 6-tuples → distance > 1.0.
        let d = ktuple_distance(b"ACDEFGHIKLMN", b"WWWWWWWWWWWW", 6);
        assert!(d > 1.0, "distance should be high for unrelated, got {d}");
    }

    #[test]
    fn ktuple_lenfac_uses_nogaplen_not_filtered_groups() {
        // C's disttbfast.c:3845-3856 computes lenfac from `nogaplen`
        // (`gappick0` strips only '-'), not from the filtered-group
        // length. For a sequence carrying 'X' / '.' / 'U', the filtered
        // group length is shorter than nogaplen; using the wrong one
        // shifts the distance by ~1e-3, which is enough to flip UPGMA
        // guide-tree join order on real-world inputs (BB12041 from
        // BALIBASE 3 is the smallest example — see `balibase_parity_run.md`).
        //
        // This test pins the behavior: two sequences with X — one short
        // enough that the filtered length differs from nogaplen — must
        // produce the same distance as the C reference computed via
        // nogaplen.
        // Two distantly related sequences so common-sextet count is
        // strictly < bunbo (i.e. the lenfac actually scales a non-zero
        // base). With X inserted, `nogaplen` increments but the filtered
        // group length is unchanged → lenfac differs → distance differs.
        let base_a = b"ACDEFGHIKLMNPQRSTVWY".repeat(8);
        let base_b = b"WYFGAILVMKRNDQECSHTP".repeat(8);

        let without_x = ktuple_distance(&base_a, &base_b, 6);

        let mut with_x_a = base_a.clone();
        let mut with_x_b = base_b.clone();
        with_x_a.insert(40, b'X');
        with_x_b.insert(40, b'X');
        let with_x = ktuple_distance(&with_x_a, &with_x_b, 6);

        // Same group content (X is filtered both times), but nogaplen
        // differs by 1 → lenfac differs → distance differs.
        let diff = (with_x - without_x).abs();
        assert!(diff > 1e-7,
            "ktuple_distance must depend on nogaplen (X kept), not filtered groups; got with_x={with_x}, without_x={without_x}, diff={diff}");
    }

    /// Sequences 69 and 649 of `many_short_3000` (lengths 29 and 30, 22
    /// shared 6-mers, self counts 24 and 25). With Apple clang's contracted
    /// `lenfac` the distance is `0x1.35450514de674p-3`; the unfused sum is
    /// one ULP lower. That ULP changed the guide tree and the alignment of
    /// the first 1726 sequences (issue #41).
    #[test]
    fn many_short_unequal_lenfac_distance_matches_c() {
        let a = b"WAWASRDQREAMMHRFTKRHTFDTMWMLV";
        let b = b"EPLWASRDQREAVMHRFTKRHTFDTMWMLV";
        let d = ktuple_distance(a, b, 6);
        let expected: u64 = if mafft_types::fp::CONTRACTS_FMA {
            0x3fc3_5450_514d_e674
        } else {
            0x3fc3_5450_514d_e673
        };
        assert_eq!(d.to_bits(), expected, "distance {:x} ({d})", d.to_bits());
    }

    /// xorshift64*, deterministic test input generator.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
        fn pick(&mut self, alphabet: &[u8]) -> u8 {
            alphabet[self.below(alphabet.len())]
        }
    }

    const PROTEIN: &[u8] = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwyXBZJUOx*.";
    const DNA: &[u8] = b"ACGTUNacgtun";

    fn random_seq(rng: &mut Rng) -> Vec<u8> {
        let len = match rng.below(10) {
            0 => rng.below(8),             // shorter than k, or empty
            1 => 6 + rng.below(40),        // short: DNA test spills into seq2
            _ => rng.below(400),
        };
        let alphabet = match rng.below(4) {
            0 => DNA,
            1 => b"ACGTACGTACGTN-" as &[u8], // nucleotide-looking with gaps
            _ => PROTEIN,
        };
        let mut s: Vec<u8> = (0..len).map(|_| rng.pick(alphabet)).collect();
        if rng.below(3) == 0 && len > 4 {
            // mostly-nucleotide prefix with one protein letter later on
            let at = rng.below(len);
            s[at] = b'W';
        }
        s
    }

    #[test]
    fn ktuple_distances_from_is_bit_identical() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for round in 0..20 {
            let seqs: Vec<Vec<u8>> = (0..40).map(|_| random_seq(&mut rng)).collect();
            for is_nucleotide in [false, true] {
                let prepared: Vec<KtupleSeq> =
                    seqs.iter().map(|s| KtupleSeq::new(s, 6, is_nucleotide)).collect();
                // one scratch reused across rows (as on a worker thread), and a
                // fresh one per row: both must give the per-pair values
                let mut shared = KtupleScratch::default();
                for i in 0..seqs.len() {
                    let row_shared = ktuple_distances_from(&prepared, i, &mut shared);
                    let row_fresh =
                        ktuple_distances_from(&prepared, i, &mut KtupleScratch::default());
                    for (off, j) in ((i + 1)..seqs.len()).enumerate() {
                        let want = ktuple_distance_typed(&seqs[i], &seqs[j], 6, is_nucleotide);
                        for got in [row_shared[off], row_fresh[off]] {
                            assert_eq!(got.to_bits(), want.to_bits(),
                                "round {round} nuc={is_nucleotide} pair ({i},{j}): {got} vs {want}\n{:?}\n{:?}",
                                String::from_utf8_lossy(&seqs[i]), String::from_utf8_lossy(&seqs[j]));
                        }
                    }
                }
            }
        }
    }

    /// The alphabet is fixed by the caller: a sequence with an ambiguity
    /// code anywhere (even among its first 100 letters, where the sniffing
    /// wrapper would have switched the pair to protein) gets the same
    /// nucleotide distance as `ktuple_distance_typed(.., true)`, and the
    /// order of the pair does not matter for the alphabet.
    #[test]
    fn ktuple_prepared_alphabet_comes_from_the_caller() {
        let seq1 = b"ACGTACGTAC".to_vec(); // 10 nucleotide letters
        for w_at in [0usize, 50, 89, 90, 95, 150] {
            let mut seq2 = b"ACGT".repeat(40);
            seq2[w_at] = b'R';
            for is_nucleotide in [true, false] {
                let fwd = [KtupleSeq::new(&seq1, 6, is_nucleotide), KtupleSeq::new(&seq2, 6, is_nucleotide)];
                let got = ktuple_distances_from(&fwd, 0, &mut KtupleScratch::default())[0];
                let want = ktuple_distance_typed(&seq1, &seq2, 6, is_nucleotide);
                assert_eq!(got.to_bits(), want.to_bits(), "nuc={is_nucleotide}, R at {w_at}");
                let rev = [KtupleSeq::new(&seq2, 6, is_nucleotide), KtupleSeq::new(&seq1, 6, is_nucleotide)];
                let got = ktuple_distances_from(&rev, 0, &mut KtupleScratch::default())[0];
                let want = ktuple_distance_typed(&seq2, &seq1, 6, is_nucleotide);
                assert_eq!(got.to_bits(), want.to_bits(), "rev, nuc={is_nucleotide}, R at {w_at}");
            }
            // The R is dropped from the nucleotide 6-mer stream wherever it
            // sits, so the prepared nucleotide sequence always has codes.
            assert!(KtupleSeq::new(&seq2, 6, true).points.is_some(), "R at {w_at}");
        }
    }

    #[test]
    fn msa_distance_prepared_is_bit_identical() {
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        // A small synthetic matrix whose amino_map also sends some letters
        // past `nalpha` (scored as "no contribution" by naive_pair_score).
        let nalpha = 6usize;
        let matrix: Vec<Vec<i32>> = (0..nalpha)
            .map(|i| (0..nalpha).map(|j| ((i * 7 + j * 3) % 11) as i32 - 4).collect())
            .collect();
        let mut amino_map = [255u8; 256];
        for (k, &c) in b"ACDEFGHIKL".iter().enumerate() {
            amino_map[c as usize] = (k % 8) as u8; // codes 6 and 7 are out of range
        }
        let real = mafft_scoring::build_context(
            mafft_types::ScoringModel::Blosum(62), mafft_types::SeqType::Protein);
        for (mtx, amap, alphabet) in [
            (&matrix, &amino_map, b"ACDEFGHIKLMNX" as &[u8]),
            (&real.substitution_matrix, &real.amino_map, PROTEIN),
        ] {
            for penalty in [-1530, 0, 7] {
                let ctx = MsaDistance::new(mtx, amap, penalty).expect("square matrix");
                for round in 0..10 {
                    let width = 1 + rng.below(300);
                    let nrow = 30;
                    let gap_rate = 1 + rng.below(9); // in tenths
                    let rows: Vec<Vec<u8>> = (0..nrow)
                        .map(|r| {
                            // a few rows of a different width exercise min(len)
                            let w = if r % 11 == 5 { rng.below(width + 20) } else { width };
                            (0..w)
                                .map(|_| if rng.below(10) < gap_rate { b'-' } else { rng.pick(alphabet) })
                                .collect()
                        })
                        .collect();
                    let prepared: Vec<ScoredRow> = rows.iter().map(|s| ctx.row(s)).collect();
                    for i in 0..nrow {
                        for j in 0..nrow {
                            let want = scoring_matrix_distance(&rows[i], &rows[j], mtx, amap, penalty);
                            let got = ctx.distance(&prepared[i], &prepared[j]);
                            assert_eq!(got.to_bits(), want.to_bits(),
                                "penalty {penalty} round {round} pair ({i},{j}): {got} vs {want}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn from_upper_rows_matches_set() {
        let mut a = DistanceMatrix::new(4);
        let mut rows = Vec::new();
        for i in 0..4 {
            let mut r = Vec::new();
            for j in (i + 1)..4 {
                let v = (i * 10 + j) as f64;
                a.set(i, j, v);
                r.push(v);
            }
            rows.push(r);
        }
        let b = DistanceMatrix::from_upper_rows(rows);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(a.get(i, j).to_bits(), b.get(i, j).to_bits());
            }
        }
    }

    #[test]
    fn half_matrix_symmetric() {
        let mut dm = DistanceMatrix::new(4);
        dm.set(0, 2, 0.5);
        assert!((dm.get(0, 2) - 0.5).abs() < 1e-10);
        assert!((dm.get(2, 0) - 0.5).abs() < 1e-10);
        assert!(dm.get(0, 0).abs() < 1e-10);
    }

    fn pseudo_dna(mut x: u64, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((x >> 33) & 3) as usize]
            })
            .collect()
    }

    /// An IUPAC code is dropped from the 6-mer stream wherever it sits
    /// (`seq_grp_nuc`, `tmp < 4`) and still counts in `nogaplen`, so the
    /// nucleotide distance must not depend on its position. The sniffing
    /// wrapper breaks this when the code is among the first 100 letters:
    /// it then takes the pair for protein (16S rRNA records with an `R`
    /// at base 68 moved the guide tree of a 250-sequence set).
    #[test]
    fn typed_nucleotide_distance_ignores_where_the_ambiguity_code_sits() {
        let other = pseudo_dna(7, 200);
        let base = pseudo_dna(11, 200);
        let mut early = base.clone();
        early.insert(10, b'R');
        let mut late = base.clone();
        late.insert(150, b'R');

        let d_early = ktuple_distance_typed(&early, &other, 6, true);
        let d_late = ktuple_distance_typed(&late, &other, 6, true);
        assert_eq!(d_early.to_bits(), d_late.to_bits(), "{d_early} vs {d_late}");
        assert_eq!(d_late.to_bits(), ktuple_distance(&late, &other, 6).to_bits());

        // The sniffing wrapper misreads `early` as protein.
        assert_ne!(ktuple_distance(&early, &other, 6).to_bits(), d_early.to_bits());
        assert_eq!(
            ktuple_distance(&early, &other, 6).to_bits(),
            ktuple_distance_typed(&early, &other, 6, false).to_bits()
        );
    }
}
