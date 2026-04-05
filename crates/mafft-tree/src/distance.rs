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
pub fn pairwise_identity_distance(seq1: &[u8], seq2: &[u8]) -> f64 {
    let mut matches = 0u64;
    let mut aligned = 0u64;

    for (a, b) in seq1.iter().zip(seq2.iter()) {
        if *a != b'-' && *b != b'-' {
            aligned += 1;
            if a == b {
                matches += 1;
            }
        }
    }

    if aligned == 0 {
        1.0
    } else {
        1.0 - matches as f64 / aligned as f64
    }
}

/// Compute k-tuple (6-mer) distance between two unaligned sequences.
///
/// Counts shared 6-tuples (sextets) between the sequences. Distance is
/// proportional to the fraction of unshared k-tuples.
///
/// This is a fast approximation used for initial guide tree construction
/// (FFT-NS-1, FFT-NS-2).
pub fn ktuple_distance(seq1: &[u8], seq2: &[u8], k: usize) -> f64 {
    if seq1.is_empty() || seq2.is_empty() || k == 0 {
        return 2.0;
    }
    if seq1.len() < k || seq2.len() < k {
        return 2.0;
    }

    // Count k-mer frequencies in seq1
    let mut table1 = std::collections::HashMap::new();
    for window in seq1.windows(k) {
        if window.iter().all(|c| c.is_ascii_alphabetic()) {
            *table1.entry(window.to_vec()).or_insert(0u32) += 1;
        }
    }

    // Count k-mer frequencies in seq2
    let mut table2 = std::collections::HashMap::new();
    for window in seq2.windows(k) {
        if window.iter().all(|c| c.is_ascii_alphabetic()) {
            *table2.entry(window.to_vec()).or_insert(0u32) += 1;
        }
    }

    // Self-scores
    let ss1: u64 = table1.values().map(|&v| v as u64).sum();
    let ss2: u64 = table2.values().map(|&v| v as u64).sum();

    if ss1 == 0 || ss2 == 0 {
        return 2.0;
    }

    // Common k-tuples: min of counts for each shared k-mer
    let common: u64 = table1
        .iter()
        .filter_map(|(kmer, &count1)| {
            table2.get(kmer).map(|&count2| count1.min(count2) as u64)
        })
        .sum();

    let min_ss = ss1.min(ss2);
    let dist = (1.0 - common as f64 / min_ss as f64) * 2.0;
    dist.clamp(0.0, 2.0)
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
        let d = ktuple_distance(b"AAAAAAAAAAAA", b"CCCCCCCCCCCC", 6);
        assert!((d - 2.0).abs() < 0.01, "distance should be ~2 for unrelated");
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
