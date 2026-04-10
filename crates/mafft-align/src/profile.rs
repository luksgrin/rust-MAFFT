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

    // Exact port of C's MSalignmm_tanni (MSalignmm.c lines 770-888).
    //
    // C's DP formulation:
    //   currentw[j] = sub(i,j) (pre-filled by match_calc)
    //   wm = max(h[i-1][j-1], mi_extend, mj_extend)
    //   h[i][j] = currentw[j] += wm  (sub added to ALL paths)
    //   mi = max(mi, h[i-1][j-1] + ogcp2[j]*gf)  (open from DIAGONAL)
    //   m[j] = max(m[j], h[i-1][j-1] + ogcp1[i]*gf)  (open from DIAGONAL)
    //
    // Traceback: ijp[i][j] = 0 (diagonal), >0 (skip rows=insertion in prof1),
    //            <0 (skip cols=deletion in prof2).

    let hgf1 = prof1.nongap_freq.first().copied().unwrap_or(1.0); // headgapfreq1
    let hgf2 = prof2.nongap_freq.first().copied().unwrap_or(1.0); // headgapfreq2
    let gf1_0 = prof1.nongap_freq.first().copied().unwrap_or(1.0); // gapfreq1f[0]
    let gf2_0 = prof2.nongap_freq.first().copied().unwrap_or(1.0); // gapfreq2f[0]

    // Full h matrix and ijp traceback (not rolling-row, but same values)
    let mut h = vec![vec![0.0f64; m + 1]; n + 1];
    let mut ijp = vec![vec![0i32; m + 1]; n + 1];

    // initverticalw[i] = match_score(i-1, 0) + head gap cost  (C lines 776, 783-786)
    // Note: C's match_calc fills initverticalw[1..n] with sub scores for
    // prof1 positions vs prof2 position 0. We compute inline.
    let mut initverticalw = vec![0.0f64; n + 1];
    for i in 1..=n {
        initverticalw[i] = prof1.match_score(i - 1, prof2, 0, matrix);
        if head_gap {
            initverticalw[i] += ogcp1[0] * hgf2 + fgcp1[i - 1] * gf2_0;
        }
    }

    // currentw[j] = match_score(0, j-1) + head gap cost  (C lines 779, 790-793)
    let mut currentw = vec![0.0f64; m + 1];
    for j in 1..=m {
        currentw[j] = prof1.match_score(0, prof2, j - 1, matrix);
        if head_gap {
            currentw[j] += ogcp2[0] * hgf1 + fgcp2[j - 1] * gf1_0;
        }
    }

    // Store h[0][j] = currentw[j] for traceback
    for j in 0..=m { h[0][j] = currentw[j]; }
    for i in 0..=n { h[i][0] = initverticalw[i]; }

    // m[j] = insertion tracker (gap in prof1), persists across rows  (C line 798)
    let mut mj = vec![f64::NEG_INFINITY; m + 1];
    let mut mpj = vec![0usize; m + 1];
    for j in 1..=m {
        let gf2_jm1 = prof2.nongap_freq.get(j - 1).copied().unwrap_or(1.0);
        mj[j] = currentw[j - 1] + ogcp1[1] * gf2_jm1;
        mpj[j] = 0;
    }

    // ijp boundary: first column = +i+1, first row = -(j+1)  (C lines 568-575)
    for i in 0..=n { ijp[i][0] = i as i32 + 1; }
    for j in 0..=m { ijp[0][j] = -(j as i32 + 1); }

    let mut previousw = vec![0.0f64; m + 1];
    let mut lastverticalw = vec![0.0f64; n + 1];
    lastverticalw[0] = currentw[m - 1];

    // Main DP loop  (C lines 807-888)
    let lasti = if tail_gap { n + 1 } else { n };
    for i in 1..lasti {
        // Swap rows  (C lines 809-811)
        std::mem::swap(&mut previousw, &mut currentw);
        previousw[0] = initverticalw[i - 1];

        // Fill currentw with match scores for row i  (C line 816)
        for j in 1..=m {
            currentw[j] = prof1.match_score(i - 1, prof2, j - 1, matrix);
        }
        currentw[0] = initverticalw[i];

        // Initialize mi (deletion tracker) for this row  (C line 819)
        let gf1_im1 = prof1.nongap_freq.get(i - 1).copied().unwrap_or(1.0);
        let mut mi = previousw[0] + ogcp2[1] * gf1_im1;
        let mut mpi: usize = 0;

        for j in 1..=m {
            let gf1_i = prof1.nongap_freq.get(i).copied().unwrap_or(1.0);
            let gf1_im1 = prof1.nongap_freq.get(i - 1).copied().unwrap_or(1.0);
            let gf2_j = prof2.nongap_freq.get(j).copied().unwrap_or(1.0);
            let gf2_jm1 = prof2.nongap_freq.get(j - 1).copied().unwrap_or(1.0);

            // wm = diagonal score  (C line 831)
            let mut wm = previousw[j - 1];
            ijp[i][j] = 0;

            // Deletion extend: mi + fgcp2[j-1] * gf1[i]  (C line 837)
            let g = mi + fgcp2[j - 1] * gf1_i;
            if g > wm {
                wm = g;
                ijp[i][j] = -(j as i32 - mpi as i32); // C line 844
            }

            // Deletion open from diagonal  (C line 846)
            let g = previousw[j - 1] + ogcp2[j] * gf1_im1;
            if g >= mi {
                mi = g;
                mpi = j - 1;
            }

            // Insertion extend: m[j] + fgcp1[i-1] * gf2[j]  (C line 856)
            let g = mj[j] + fgcp1[i - 1] * gf2_j;
            if g > wm {
                wm = g;
                ijp[i][j] = i as i32 - mpj[j] as i32; // C line 863
            }

            // Insertion open from diagonal  (C line 865)
            let g = previousw[j - 1] + ogcp1[i] * gf2_jm1;
            if g >= mj[j] {
                mj[j] = g;
                mpj[j] = i - 1;
            }

            // h[i][j] = sub(i,j) + wm  (C line 878: *curpt += wm)
            currentw[j] += wm;
            h[i][j] = currentw[j];
        }
        lastverticalw[i] = currentw[m - 1];
    }

    // Tail gap handling  (C lines 512-536)
    if !tail_gap {
        let mut wm = lastverticalw[0];
        for i in 0..n {
            if lastverticalw[i] >= wm {
                wm = lastverticalw[i];
                ijp[n][m] = (n - i) as i32;
            }
        }
        for j in 0..m {
            if h[n - 1][j] >= wm {
                wm = h[n - 1][j];
                ijp[n][m] = -((m - j) as i32);
            }
        }
    }

    // Traceback  (C lines 581-617)
    let mut gaptable1 = Vec::new(); // 'o' = content, '-' = gap
    let mut gaptable2 = Vec::new();

    let (mut iin, mut jin) = (n as i32, m as i32);
    let klim = n + m;
    let mut k = 0;
    while k <= klim {
        let (ifi, jfi): (i32, i32);
        let v = ijp[iin as usize][jin as usize];
        if v < 0 {
            // Deletion: skip columns
            ifi = iin - 1;
            jfi = jin + v;
        } else if v > 0 {
            // Insertion: skip rows
            ifi = iin - v;
            jfi = jin - 1;
        } else {
            // Diagonal
            ifi = iin - 1;
            jfi = jin - 1;
        }

        // Emit gap in prof2 for skipped rows (insertion in prof1)
        let mut l = iin - ifi;
        while l > 1 {
            gaptable1.push(b'o');
            gaptable2.push(b'-');
            k += 1;
            l -= 1;
        }
        // Emit gap in prof1 for skipped columns (deletion in prof2)
        let mut l = jin - jfi;
        while l > 1 {
            gaptable1.push(b'-');
            gaptable2.push(b'o');
            k += 1;
            l -= 1;
        }

        if iin <= 0 || jin <= 0 { break; }
        // Emit diagonal match
        gaptable1.push(b'o');
        gaptable2.push(b'o');
        k += 1;
        iin = ifi;
        jin = jfi;
    }

    // Convert gaptable (built in reverse) to AlignOps
    gaptable1.reverse();
    gaptable2.reverse();

    let mut ops = Vec::with_capacity(gaptable1.len());
    for k in 0..gaptable1.len() {
        match (gaptable1[k], gaptable2[k]) {
            (b'o', b'o') => ops.push(AlignOp::Match),
            (b'o', b'-') => ops.push(AlignOp::Delete),  // gap in prof2
            (b'-', b'o') => ops.push(AlignOp::Insert),  // gap in prof1
            _ => {}
        }
    }

    let best_score = h[n][m];

    Alignment {
        seq1: Vec::new(),
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
