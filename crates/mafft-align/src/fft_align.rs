/// FFT-accelerated alignment.
///
/// Ports the C `Falign()` from Falign.c.
///
/// Uses FFT cross-correlation to find anchor points between two sequence
/// groups, then applies full DP alignment within each anchored segment.

use mafft_fft::{cross_correlate, get_top_candidates, alignable_segments, SegmentParams};

use crate::dp::{Alignment, AlignOp, GapModel};
use crate::profile::{Profile, profile_align};

/// Parameters controlling FFT-accelerated alignment.
#[derive(Debug, Clone)]
pub struct FftAlignParams {
    /// Number of candidate lags to evaluate.
    pub num_candidates: usize,
    /// Segment detection parameters.
    pub segment_params: SegmentParams,
    /// Gap model for DP alignment within segments.
    pub gap: GapModel,
    /// Whether to penalize head/tail gaps.
    pub head_gap: bool,
    pub tail_gap: bool,
}

impl FftAlignParams {
    pub fn protein() -> Self {
        Self {
            num_candidates: 20,
            segment_params: SegmentParams::protein(),
            gap: GapModel::default(),
            head_gap: true,
            tail_gap: true,
        }
    }

    pub fn dna() -> Self {
        Self {
            num_candidates: 20,
            segment_params: SegmentParams::dna(),
            gap: GapModel::default(),
            head_gap: true,
            tail_gap: true,
        }
    }
}

/// An anchor point between two profiles.
#[derive(Debug, Clone)]
pub struct Anchor {
    /// Position in profile 1.
    pub pos1: usize,
    /// Position in profile 2.
    pub pos2: usize,
}

/// Perform FFT-accelerated profile alignment.
///
/// 1. Converts profiles to numerical vectors.
/// 2. Uses FFT cross-correlation to find the best lag (shift).
/// 3. Detects alignable segments at that lag.
/// 4. Runs full DP alignment within each segment.
/// 5. Concatenates segment alignments.
///
/// Falls back to direct profile DP if no good anchors are found.
pub fn fft_profile_align(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    params: &FftAlignParams,
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

    // Step 1: Build score vectors for correlation.
    // For each position, compute the diagonal score if the two profiles
    // were aligned at that offset. We sum across the alphabet.
    let site_scores = compute_site_scores(prof1, prof2, matrix);

    // Step 2: Find best lag via FFT correlation
    let vec1 = profile_to_score_vector(prof1, matrix);
    let vec2 = profile_to_score_vector(prof2, matrix);
    let correlation = cross_correlate(&vec1, &vec2);
    let candidates = get_top_candidates(&correlation, params.num_candidates);

    // Step 3: For the best candidate, detect alignable segments
    let mut best_segments: Vec<mafft_fft::AlignableSegment> = Vec::new();
    let mut best_lag = 0i32;

    for cand in &candidates {
        let lag = cand.lag;
        let shifted_scores = shift_and_score(prof1, prof2, matrix, lag);
        let segments = alignable_segments(&shifted_scores, &params.segment_params);
        if segments.len() > best_segments.len()
            || (segments.len() == best_segments.len()
                && segments.iter().map(|s| s.score).sum::<f64>()
                    > best_segments.iter().map(|s| s.score).sum::<f64>())
        {
            best_segments = segments;
            best_lag = lag;
        }
    }

    // Step 4: If no anchors found, fall back to direct DP
    if best_segments.is_empty() {
        return profile_align(prof1, prof2, matrix, &params.gap, params.head_gap, params.tail_gap);
    }

    // Step 5: Convert segments to anchors and align between them
    let anchors: Vec<Anchor> = best_segments
        .iter()
        .map(|seg| {
            let center = seg.center;
            if best_lag >= 0 {
                Anchor {
                    pos1: center,
                    pos2: (center as i32 - best_lag).max(0) as usize,
                }
            } else {
                Anchor {
                    pos1: (center as i32 + best_lag).max(0) as usize,
                    pos2: center,
                }
            }
        })
        .filter(|a| a.pos1 < n && a.pos2 < m)
        .collect();

    if anchors.is_empty() {
        return profile_align(prof1, prof2, matrix, &params.gap, params.head_gap, params.tail_gap);
    }

    // Step 6: Align segments between anchors using full DP
    align_with_anchors(prof1, prof2, matrix, &params.gap, &anchors)
}

/// Align two profiles using pre-computed anchor points.
///
/// Divides the alignment into segments defined by anchors and aligns
/// each segment independently, then concatenates.
fn align_with_anchors(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    gap: &GapModel,
    anchors: &[Anchor],
) -> Alignment {
    let mut all_ops = Vec::new();
    let mut total_score = 0.0;

    let mut p1 = 0usize;
    let mut p2 = 0usize;

    for anchor in anchors {
        // Align the segment before this anchor
        if anchor.pos1 > p1 || anchor.pos2 > p2 {
            let sub1 = sub_profile(prof1, p1, anchor.pos1);
            let sub2 = sub_profile(prof2, p2, anchor.pos2);
            let seg_aln = profile_align(&sub1, &sub2, matrix, gap, p1 == 0, false);
            total_score += seg_aln.score;
            all_ops.extend(seg_aln.operations);
        }

        // The anchor itself is a match
        if anchor.pos1 < prof1.length && anchor.pos2 < prof2.length {
            all_ops.push(AlignOp::Match);
            p1 = anchor.pos1 + 1;
            p2 = anchor.pos2 + 1;
        }
    }

    // Align the tail after the last anchor
    if p1 < prof1.length || p2 < prof2.length {
        let sub1 = sub_profile(prof1, p1, prof1.length);
        let sub2 = sub_profile(prof2, p2, prof2.length);
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

/// Extract a sub-profile (slice of positions).
fn sub_profile(prof: &Profile, start: usize, end: usize) -> Profile {
    let end = end.min(prof.length);
    let start = start.min(end);
    Profile {
        freqs: prof.freqs[start..end].to_vec(),
        gap_freq: prof.gap_freq[start..end].to_vec(),
        length: end - start,
        nalphabets: prof.nalphabets,
    }
}

/// Convert a profile to a single score vector for FFT correlation.
///
/// Uses the diagonal self-score at each position (sum of freq² × matrix diagonal).
fn profile_to_score_vector(prof: &Profile, matrix: &[Vec<i32>]) -> Vec<f64> {
    (0..prof.length)
        .map(|pos| {
            let mut s = 0.0;
            for a in 0..prof.nalphabets.min(matrix.len()) {
                s += prof.freqs[pos][a] * matrix[a][a] as f64 * prof.freqs[pos][a];
            }
            s
        })
        .collect()
}

/// Compute per-position match scores between two profiles at a given lag.
fn shift_and_score(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    lag: i32,
) -> Vec<f64> {
    let len = prof1.length.max(prof2.length);
    let mut scores = vec![0.0; len];

    for i in 0..prof1.length {
        let j = i as i32 - lag;
        if j >= 0 && (j as usize) < prof2.length {
            scores[i] = prof1.match_score(i, prof2, j as usize, matrix);
        }
    }
    scores
}

/// Compute site scores between two profiles without lag (for direct use).
fn compute_site_scores(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
) -> Vec<f64> {
    shift_and_score(prof1, prof2, matrix, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_setup() -> (Vec<Vec<i32>>, [u8; 256], usize) {
        let mut mtx = vec![vec![-100i32; 5]; 5];
        for i in 0..4 { mtx[i][i] = 100; }
        let mut map = [0xFFu8; 256];
        map[b'A' as usize] = 0;
        map[b'C' as usize] = 1;
        map[b'G' as usize] = 2;
        map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;
        (mtx, map, 5)
    }

    #[test]
    fn fft_align_identical_profiles() {
        let (mtx, map, nalpha) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGTACGTACGTACGTACGTACGT"];
        let prof = Profile::from_aligned(&seqs, &[1.0], &map, nalpha);
        let params = FftAlignParams {
            num_candidates: 5,
            segment_params: SegmentParams {
                window_size: 4,
                threshold: 100.0,
                max_segment_size: 50,
            },
            gap: GapModel::new(-200.0, -10.0),
            head_gap: true,
            tail_gap: true,
        };
        let aln = fft_profile_align(&prof, &prof, &mtx, &params);
        // Should produce a good alignment (exact score depends on anchor placement)
        assert!(aln.score > 0.0);
        assert!(!aln.operations.is_empty());
    }

    #[test]
    fn fft_align_falls_back_on_short() {
        let (mtx, map, nalpha) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGT"];
        let prof = Profile::from_aligned(&seqs, &[1.0], &map, nalpha);
        let params = FftAlignParams::protein();
        // Short sequences should fall back to direct DP
        let aln = fft_profile_align(&prof, &prof, &mtx, &params);
        assert!((aln.score - 400.0).abs() < 1e-6);
    }
}
