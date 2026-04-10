/// Profile-to-profile (group-to-group) alignment.
///
/// Ports the C `MSalignmm()` from MSalignmm.c.
///
/// Aligns two groups of sequences by computing position-specific frequency
/// matrices (profiles) and running affine-gap DP on the profile scores.

use crate::dp::{AlignOp, Alignment, GapModel};

/// A position-specific frequency matrix (profile).
///
/// `freq[position][alphabet_index]` = weighted frequency of residue at that position.
#[derive(Debug, Clone)]
pub struct Profile {
    /// Frequency matrix: `freqs[pos][residue_idx]`.
    pub freqs: Vec<Vec<f64>>,
    /// Gap frequency at each position (0.0 = no gaps, 1.0 = all gaps).
    pub gap_freq: Vec<f64>,
    /// Non-gap frequency at each position (= 1.0 - gap_freq). Used for
    /// cross-weighting in gap cost computation.
    pub nongap_freq: Vec<f64>,
    /// Opening gap cost profile: weighted count of gap-opening transitions
    /// (non-gap → gap) at each position, modulated by penalty and nongap_freq.
    /// Formula: `ogcp[i] = 0.5 * (1.0 - opening_count[i]) * penalty * nongap_freq[i]`
    pub ogcp: Vec<f64>,
    /// Final gap cost profile: weighted count of gap-closing transitions
    /// (gap → non-gap) at each position.
    /// Formula: `fgcp[i] = 0.5 * (1.0 - closing_count[i]) * penalty * nongap_freq[i]`
    pub fgcp: Vec<f64>,
    /// Number of positions (alignment columns).
    pub length: usize,
    /// Alphabet size.
    pub nalphabets: usize,
}

impl Profile {
    /// Build a profile from a group of aligned sequences with weights.
    ///
    /// - `sequences`: aligned sequences (same length, with gaps as '-').
    /// - `weights`: per-sequence weights (must sum to ~1.0 for best results).
    /// - `amino_map`: ASCII char → internal residue index.
    /// - `nalphabets`: alphabet size (20 for protein, 26 with ambiguity).
    pub fn from_aligned(
        sequences: &[&[u8]],
        weights: &[f64],
        amino_map: &[u8; 256],
        nalphabets: usize,
    ) -> Self {
        assert_eq!(sequences.len(), weights.len());
        if sequences.is_empty() {
            return Self {
                freqs: Vec::new(),
                gap_freq: Vec::new(),
                nongap_freq: Vec::new(),
                ogcp: Vec::new(),
                fgcp: Vec::new(),
                length: 0,
                nalphabets,
            };
        }

        let length = sequences[0].len();
        let mut freqs = vec![vec![0.0f64; nalphabets]; length];
        let mut gap_freq = vec![0.0f64; length];

        // Opening gap count (non-gap → gap transitions) and final gap count
        // (gap → non-gap transitions), matching C's st_OpeningGapCount/st_FinalGapCount.
        let mut opening_count = vec![0.0f64; length];
        let mut closing_count = vec![0.0f64; length];

        for (seq, &w) in sequences.iter().zip(weights.iter()) {
            let mut prev_is_gap = false;
            for (pos, &ch) in seq.iter().enumerate() {
                if pos >= length { break; }
                let is_gap = ch == b'-' || ch == b'.';
                if is_gap {
                    gap_freq[pos] += w;
                } else {
                    let idx = amino_map[ch as usize] as usize;
                    if idx < nalphabets {
                        freqs[pos][idx] += w;
                    }
                }
                // Gap opening: non-gap → gap
                if !prev_is_gap && is_gap {
                    opening_count[pos] += w;
                }
                // Gap closing: gap → non-gap
                if prev_is_gap && !is_gap {
                    closing_count[pos] += w;
                }
                prev_is_gap = is_gap;
            }
        }

        let nongap_freq: Vec<f64> = gap_freq.iter().map(|&g| (1.0 - g).max(0.0)).collect();

        // Clamp opening/closing counts to [0, 1] range.
        // C's st_OpeningGapCount/st_FinalGapCount produce values in this
        // range because sequence weights sum to 1. Our weights may not be
        // normalized, so we clamp.
        for v in &mut opening_count { *v = v.clamp(0.0, 1.0); }
        for v in &mut closing_count { *v = v.clamp(0.0, 1.0); }

        // Store raw opening/closing counts; actual ogcp/fgcp are computed
        // in profile_align when the penalty parameter is known.
        Self {
            freqs,
            gap_freq,
            nongap_freq,
            ogcp: opening_count,
            fgcp: closing_count,
            length,
            nalphabets,
        }
    }

    /// Compute the match score between position `i` of this profile and
    /// position `j` of another profile, using the given scoring matrix.
    ///
    /// Uses a two-pass approach (matching C's `match_calc`) that enables
    /// SIMD auto-vectorization:
    /// 1. Build `scarr[b] = sum_a(freq1[a] * matrix[a][b])` — branchless
    /// 2. Dot-product `sum_b(scarr[b] * freq2[b])` — branchless, contiguous
    #[inline]
    pub fn match_score(
        &self,
        i: usize,
        other: &Profile,
        j: usize,
        matrix: &[Vec<i32>],
    ) -> f64 {
        let nalpha = self.nalphabets.min(other.nalphabets).min(matrix.len());
        let freq1 = &self.freqs[i];
        let freq2 = &other.freqs[j];

        // Pass 1: scarr[b] = sum_a(freq1[a] * matrix[a][b])
        // This is branchless — zero freq1 values contribute zero, no skip needed.
        let mut scarr = [0.0f64; 32]; // fixed-size for auto-vectorization (covers nalphabets <= 26)
        for a in 0..nalpha {
            let f1 = freq1[a];
            // Compiler can vectorize this inner loop: scarr[b] += f1 * matrix[a][b]
            let row = &matrix[a];
            let row_len = nalpha.min(row.len());
            for b in 0..row_len {
                scarr[b] += f1 * row[b] as f64;
            }
        }

        // Pass 2: dot product — perfectly vectorizable contiguous f64 multiply-add
        let mut score = 0.0f64;
        for b in 0..nalpha {
            score += scarr[b] * freq2[b];
        }
        score
    }

    /// Extract a sub-profile (slice of positions from `start` to `end`).
    pub fn sub_profile(&self, start: usize, end: usize) -> Profile {
        let end = end.min(self.length);
        let start = start.min(end);
        Profile {
            freqs: self.freqs[start..end].to_vec(),
            gap_freq: self.gap_freq[start..end].to_vec(),
            nongap_freq: self.nongap_freq[start..end].to_vec(),
            ogcp: self.ogcp[start..end].to_vec(),
            fgcp: self.fgcp[start..end].to_vec(),
            length: end - start,
            nalphabets: self.nalphabets,
        }
    }
}

/// Align two profiles using anchor points, running DP within each segment.
///
/// Shared implementation used by both `fft_align` and `constrained_align`.
pub fn align_with_anchors(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    gap: &GapModel,
    anchors: &[(usize, usize)],
) -> Alignment {
    let mut all_ops = Vec::new();
    let mut total_score = 0.0;
    let mut p1 = 0usize;
    let mut p2 = 0usize;

    for &(a1, a2) in anchors {
        if a1 > p1 || a2 > p2 {
            let sub1 = prof1.sub_profile(p1, a1);
            let sub2 = prof2.sub_profile(p2, a2);
            let seg_aln = profile_align(&sub1, &sub2, matrix, gap, p1 == 0, false);
            total_score += seg_aln.score;
            all_ops.extend(seg_aln.operations);
        }
        if a1 < prof1.length && a2 < prof2.length {
            all_ops.push(AlignOp::Match);
            p1 = a1 + 1;
            p2 = a2 + 1;
        }
    }

    if p1 < prof1.length || p2 < prof2.length {
        let sub1 = prof1.sub_profile(p1, prof1.length);
        let sub2 = prof2.sub_profile(p2, prof2.length);
        let seg_aln = profile_align(&sub1, &sub2, matrix, gap, false, true);
        total_score += seg_aln.score;
        all_ops.extend(seg_aln.operations);
    }

    Alignment {
        seq1: Vec::new(),
        seq2: Vec::new(),
        score: total_score,
        operations: all_ops,
    }
}

/// Align two profiles using affine gap DP.
///
/// Gap penalties are modulated by gap frequency: positions where many
/// sequences already have gaps receive reduced gap opening penalties.
///
/// Returns the alignment operations and score.
pub fn profile_align(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    gap: &GapModel,
    head_gap: bool,
    tail_gap: bool,
) -> Alignment {
    let n = prof1.length;
    let m = prof2.length;

    if n == 0 || m == 0 {
        return Alignment {
            seq1: Vec::new(),
            seq2: Vec::new(),
            score: 0.0,
            operations: Vec::new(),
        };
    }

    // Compute position-specific gap cost profiles matching C's formula:
    // ogcp[i] = 0.5 * (1.0 - opening_count[i]) * penalty * nongap_freq[i]
    // fgcp[i] = 0.5 * (1.0 - closing_count[i]) * penalty * nongap_freq[i]
    //
    // C allocates these arrays with size lgth+2 (via AllocateFloatVec),
    // so positions beyond the profile length are 0.0 (from calloc).
    // The DP loop accesses ogcp2[j] where j goes 1..lgth2, so the last
    // element accessed is ogcp2[lgth2] = 0.0. We add a trailing 0.0
    // to match this behavior.
    let penalty = gap.open;
    let mut ogcp1: Vec<f64> = (0..n).map(|i| {
        0.5 * (1.0 - prof1.ogcp[i]) * penalty * prof1.nongap_freq[i]
    }).collect();
    ogcp1.push(0.0); // C's calloc padding
    let mut fgcp1: Vec<f64> = (0..n).map(|i| {
        0.5 * (1.0 - prof1.fgcp[i]) * penalty * prof1.nongap_freq[i]
    }).collect();
    fgcp1.push(0.0);
    let mut ogcp2: Vec<f64> = (0..m).map(|j| {
        0.5 * (1.0 - prof2.ogcp[j]) * penalty * prof2.nongap_freq[j]
    }).collect();
    ogcp2.push(0.0);
    let mut fgcp2: Vec<f64> = (0..m).map(|j| {
        0.5 * (1.0 - prof2.fgcp[j]) * penalty * prof2.nongap_freq[j]
    }).collect();
    fgcp2.push(0.0);

    let head_factor = if head_gap { 1.0 } else { 0.0 };

    // DP matrices
    let mut h = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut d = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut ins = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut traceback = vec![vec![0u8; m + 1]; n + 1];

    h[0][0] = 0.0;
    for i in 1..=n {
        // Terminal gap: ogcp1[0] * headgapfreq2 + fgcp1[i-1] * nongap2[0]
        let cost = if head_gap {
            ogcp1[0] * prof2.nongap_freq.first().copied().unwrap_or(1.0)
                + fgcp1[i - 1] * prof2.nongap_freq.first().copied().unwrap_or(1.0)
        } else {
            gap.open * 0.0 // TERMGAPFAC = 0 for free terminal gaps
        };
        d[i][0] = if i == 1 { cost } else { d[i-1][0] + fgcp1[i-1] * prof2.nongap_freq.first().copied().unwrap_or(1.0) };
        h[i][0] = d[i][0];
        traceback[i][0] = 1;
    }
    for j in 1..=m {
        let cost = if head_gap {
            ogcp2[0] * prof1.nongap_freq.first().copied().unwrap_or(1.0)
                + fgcp2[j - 1] * prof1.nongap_freq.first().copied().unwrap_or(1.0)
        } else {
            0.0
        };
        ins[0][j] = if j == 1 { cost } else { ins[0][j-1] + fgcp2[j-1] * prof1.nongap_freq.first().copied().unwrap_or(1.0) };
        h[0][j] = ins[0][j];
        traceback[0][j] = 2;
    }

    for i in 1..=n {
        for j in 1..=m {
            let sub = prof1.match_score(i - 1, prof2, j - 1, matrix);

            let diag = h[i - 1][j - 1] + sub;

            // Gap costs using C's exact formula (MSalignmm.c lines 837-865):
            // Deletion (gap in seq2):
            //   extend: fgcp2[j-1] * gapfreq1f[i]   (C line 837)
            //   open:   ogcp2[j]   * gapfreq1f[i-1]  (C line 846)
            // Insertion (gap in seq1):
            //   extend: fgcp1[i-1] * gapfreq2f[j]   (C line 856)
            //   open:   ogcp1[i]   * gapfreq2f[j-1]  (C line 865)
            //
            // In C, i,j are 1-based. gapfreq = nongap_freq (0-based).
            // ogcp/fgcp arrays have length lgth+2 (indices 0..lgth+1).
            // C's gapfreq arrays have size lgth+1, computed via gapcountf then
            // inverted (1.0 - gap_count). The element at index lgth = 1.0
            // (no gap at the tail position). Use 1.0 as out-of-bounds fallback.
            let gf1_i = prof1.nongap_freq.get(i).copied().unwrap_or(1.0);     // gapfreq1f[i]
            let gf1_im1 = prof1.nongap_freq.get(i - 1).copied().unwrap_or(1.0); // gapfreq1f[i-1]
            let gf2_j = prof2.nongap_freq.get(j).copied().unwrap_or(1.0);     // gapfreq2f[j]
            let gf2_jm1 = prof2.nongap_freq.get(j - 1).copied().unwrap_or(1.0); // gapfreq2f[j-1]

            let d_ext = d[i - 1][j] + fgcp2[j - 1] * gf1_i;
            let d_open = h[i - 1][j] + ogcp2[j] * gf1_im1;
            d[i][j] = d_ext.max(d_open);

            let i_ext = ins[i][j - 1] + fgcp1[i - 1] * gf2_j;
            let i_open = h[i][j - 1] + ogcp1[i] * gf2_jm1;
            ins[i][j] = i_ext.max(i_open);

            h[i][j] = diag;
            traceback[i][j] = 0;
            if d[i][j] > h[i][j] { h[i][j] = d[i][j]; traceback[i][j] = 1; }
            if ins[i][j] > h[i][j] { h[i][j] = ins[i][j]; traceback[i][j] = 2; }
        }
    }

    // Find best endpoint
    let (mut ei, mut ej) = (n, m);
    let mut best_score = h[n][m];

    if !tail_gap {
        for j in 0..m {
            if h[n][j] > best_score { best_score = h[n][j]; ei = n; ej = j; }
        }
        for i in 0..n {
            if h[i][m] > best_score { best_score = h[i][m]; ei = i; ej = m; }
        }
    }

    // Traceback
    let mut ops = Vec::new();
    let (mut i, mut j) = (ei, ej);

    while i < n { ops.push(AlignOp::Delete); i += 1; }
    while j < m { ops.push(AlignOp::Insert); j += 1; }

    i = ei;
    j = ej;
    while i > 0 || j > 0 {
        match traceback[i][j] {
            0 if i > 0 && j > 0 => { ops.push(AlignOp::Match); i -= 1; j -= 1; }
            1 if i > 0 => { ops.push(AlignOp::Delete); i -= 1; }
            2 if j > 0 => { ops.push(AlignOp::Insert); j -= 1; }
            _ => break,
        }
    }
    ops.reverse();

    Alignment {
        seq1: Vec::new(), // profiles don't produce sequence strings directly
        seq2: Vec::new(),
        score: best_score,
        operations: ops,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_setup() -> (Vec<Vec<i32>>, [u8; 256]) {
        let mut mtx = vec![vec![-100i32; 5]; 5];
        for i in 0..4 { mtx[i][i] = 100; }
        let mut map = [0xFFu8; 256];
        map[b'A' as usize] = 0;
        map[b'C' as usize] = 1;
        map[b'G' as usize] = 2;
        map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;
        (mtx, map)
    }

    #[test]
    fn profile_from_single_sequence() {
        let (_, map) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGT"];
        let prof = Profile::from_aligned(&seqs, &[1.0], &map, 5);
        assert_eq!(prof.length, 4);
        assert!((prof.freqs[0][0] - 1.0).abs() < 1e-10); // A at pos 0
        assert!((prof.freqs[1][1] - 1.0).abs() < 1e-10); // C at pos 1
    }

    #[test]
    fn profile_match_score_identical() {
        let (mtx, map) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGT"];
        let prof = Profile::from_aligned(&seqs, &[1.0], &map, 5);
        // Score of position 0 vs position 0 should be 100 (A vs A)
        let score = prof.match_score(0, &prof, 0, &mtx);
        assert!((score - 100.0).abs() < 1e-6);
    }

    #[test]
    fn profile_align_identical() {
        let (mtx, map) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGT"];
        let prof = Profile::from_aligned(&seqs, &[1.0], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);
        let aln = profile_align(&prof, &prof, &mtx, &gap, true, true);
        assert!((aln.score - 400.0).abs() < 1e-6);
        assert_eq!(aln.operations.len(), 4);
        assert!(aln.operations.iter().all(|op| *op == AlignOp::Match));
    }

    #[test]
    fn profile_with_gaps_reduces_penalty() {
        let (_mtx, map) = simple_setup();
        // Two sequences, one has a gap at position 1
        let seqs: Vec<&[u8]> = vec![b"A-GT", b"ACGT"];
        let prof = Profile::from_aligned(&seqs, &[0.5, 0.5], &map, 5);
        assert!((prof.gap_freq[1] - 0.5).abs() < 1e-10); // 50% gap at pos 1
    }
}
