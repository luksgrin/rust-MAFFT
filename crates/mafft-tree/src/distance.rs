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

/// Compute k-tuple (6-mer) distance between two unaligned sequences.
///
/// Ports the C `commonsextet_p` + `distcompact` algorithm exactly:
/// - Amino acids are mapped to 6 physicochemical groups (matching C's `amino_grp`)
/// - DNA bases are mapped to 4 groups
/// - K-tuples are encoded as integers into a flat frequency table
/// - Distance includes C's length adjustment factor (`lenfac`)
///
/// Formula: `(1.0 - common/min(ss1,ss2)) * lenfac * 2.0`
pub fn ktuple_distance(seq1: &[u8], seq2: &[u8], k: usize) -> f64 {
    if seq1.is_empty() || seq2.is_empty() || k == 0 {
        return 2.0;
    }
    if seq1.len() < k || seq2.len() < k {
        return 2.0;
    }

    // Detect if DNA or protein based on content
    let is_dna = seq1.iter().chain(seq2.iter())
        .filter(|&&c| c.is_ascii_alphabetic())
        .take(100)
        .all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u' | b'N' | b'n'));

    if is_dna {
        ktuple_distance_nuc(seq1, seq2, k)
    } else {
        ktuple_distance_aa(seq1, seq2, k)
    }
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

    // Length factor (C's protein parameters)
    let len1 = groups1.len() as f64;
    let len2 = groups2.len() as f64;
    let lenfac = compute_lenfac(len1, len2, 0.01, 10000.0, 10000.0, 0.1);

    let dist = (1.0 - common as f64 / ss1.min(ss2) as f64) * lenfac * 2.0;
    dist.clamp(0.0, 2.0)
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

    // Length factor (C's DNA 6-tuple parameters)
    let len1 = groups1.len() as f64;
    let len2 = groups2.len() as f64;
    let lenfac = compute_lenfac(len1, len2, 0.01, 2500.0, 2500.0, 0.1);

    let dist = (1.0 - common as f64 / ss1.min(ss2) as f64) * lenfac * 2.0;
    dist.clamp(0.0, 2.0)
}

/// C's length adjustment factor.
///
/// `lenfac = 1.0 / (shorter/longer * d + b/(longer + c) + a)`
fn compute_lenfac(len1: f64, len2: f64, a: f64, b: f64, c: f64, d: f64) -> f64 {
    let longer = len1.max(len2);
    let shorter = len1.min(len2);
    if longer == 0.0 {
        return 1.0;
    }
    1.0 / (shorter / longer * d + b / (longer + c) + a)
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
        // Use clearly protein sequences (W, F are not DNA bases)
        let d = ktuple_distance(b"WWWWWWWWWWWW", b"FFFFFFFFFFF", 6);
        // Both map to aromatic group 4, but different amino acids → same group
        // So use sequences from different groups:
        let d = ktuple_distance(b"ACDEFGHIKLMN", b"WWWWWWWWWWWW", 6);
        // No common 6-tuples (mixed groups vs all-aromatic) → distance > 1.0
        assert!(d > 1.0, "distance should be high for unrelated, got {d}");
    }

    #[test]
    fn half_matrix_symmetric() {
        let mut dm = DistanceMatrix::new(4);
        dm.set(0, 2, 0.5);
        assert!((dm.get(0, 2) - 0.5).abs() < 1e-10);
        assert!((dm.get(2, 0) - 0.5).abs() < 1e-10);
        assert!(dm.get(0, 0).abs() < 1e-10);
    }
}
