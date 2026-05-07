/// FFT-accelerated alignment.
///
/// Ports the C `Falign()` from Falign.c.
///
/// Uses multi-channel FFT cross-correlation to find anchor points between
/// two sequence groups, selects optimal anchors via DP, then applies full
/// DP alignment within each anchored segment.

use mafft_fft::{
    alignable_segments, block_align, get_top_candidates,
    multichannel_correlate, AlignableSegment, SegmentParams,
};

use crate::dp::{Alignment, GapModel};
use crate::profile::{align_with_anchors, profile_align, Profile};

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
    /// Number of FFT channels (20 for protein, 4 for DNA).
    pub num_channels: usize,
}

impl FftAlignParams {
    pub fn protein() -> Self {
        Self {
            num_candidates: 20,
            segment_params: SegmentParams::protein(),
            gap: GapModel::default(),
            head_gap: true,
            tail_gap: true,
            num_channels: 20,
        }
    }

    pub fn dna() -> Self {
        Self {
            num_candidates: 20,
            segment_params: SegmentParams::dna(),
            gap: GapModel::default(),
            head_gap: true,
            tail_gap: true,
            num_channels: 4,
        }
    }
}

/// An anchor point between two profiles.
#[derive(Debug, Clone)]
pub struct Anchor {
    pub pos1: usize,
    pub pos2: usize,
}

/// Find FFT-based anchor points between two profiles.
///
/// Performs the FFT correlation, segment detection, and anchor selection
/// steps of Falign, returning the anchors without running the DP.
/// Returns `None` if no good anchors are found.
pub fn find_fft_anchors(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    params: &FftAlignParams,
) -> Option<Vec<(usize, usize)>> {
    let n = prof1.length;
    let m = prof2.length;

    if n == 0 || m == 0 {
        return None;
    }

    let fft_size = n.max(m).next_power_of_two();
    let channels_a = profile_to_channels(prof1, params.num_channels, fft_size);
    let channels_b = profile_to_channels(prof2, params.num_channels, fft_size);
    let raw_corr = multichannel_correlate(&channels_a, &channels_b);

    let nlen = fft_size;
    let nlen2 = nlen / 2;
    let mut soukan = vec![0.0f64; nlen];
    for idx in 0..nlen {
        let lag = idx as i32 - nlen2 as i32;
        let src = lag.rem_euclid(nlen as i32) as usize;
        soukan[idx] = raw_corr[src];
    }

    let candidates = get_top_candidates(&soukan, params.num_candidates);

    // Mirrors C's `Falign` (`Falign.c:1307-1372`): accumulate segments from
    // candidate lags into a single flat list (sorted by score), then
    // independently sort centers in each seq dimension and run
    // `blockAlign2` over a sparse cross-score matrix. The DP picks the
    // best non-conflicting subset of anchors across all candidate lags.
    //
    // C's loop bounds: `if (lag <= -len1 || lag >= len2) continue;`
    // (`Falign.c:1310`) and break on first empty `tmpint == 0`
    // (`Falign.c:1330`).
    struct PairSeg { center1: usize, center2: usize, score: f64 }
    let mut all: Vec<PairSeg> = Vec::new();
    for cand in &candidates {
        let lag = cand.lag;
        if lag <= -(n as i32) || lag >= m as i32 { continue; }
        let shifted_scores = shift_and_score(prof1, prof2, matrix, lag);
        let segments = alignable_segments(&shifted_scores, &params.segment_params);
        if segments.is_empty() {
            break; // C: `if(tmpint == 0) break;`
        }
        for seg in segments {
            // shift_and_score puts scores[i] = match(prof1[i], prof2[i - lag]).
            // Segment.center maps to (seq1=center, seq2=center-lag) — the
            // sign is opposite to C's because our FFT correlation peak
            // convention is mirrored relative to C's `zurasu2`. Empirically
            // this matches all FFT-NS-2 byte tests.
            let c1 = seg.center;
            let c2_signed = seg.center as i32 - lag;
            if c1 >= n || c2_signed < 0 || (c2_signed as usize) >= m { continue; }
            all.push(PairSeg { center1: c1, center2: c2_signed as usize, score: seg.score });
        }
    }

    if all.is_empty() { return None; }

    let nseg = all.len();
    let mut sort1: Vec<usize> = (0..nseg).collect();
    sort1.sort_by(|&a, &b| all[a].center1.cmp(&all[b].center1));
    let mut sort2: Vec<usize> = (0..nseg).collect();
    sort2.sort_by(|&a, &b| all[a].center2.cmp(&all[b].center2));

    let mut rank1 = vec![0usize; nseg];
    let mut rank2 = vec![0usize; nseg];
    for (r, &i) in sort1.iter().enumerate() { rank1[i] = r; }
    for (r, &i) in sort2.iter().enumerate() { rank2[i] = r; }

    let size = nseg + 2;
    let mut crossscore = vec![vec![0.0f64; size]; size];
    for i in 0..nseg {
        crossscore[rank1[i] + 1][rank2[i] + 1] = all[i].score;
    }
    crossscore[0][0] = 1e7;
    crossscore[size - 1][size - 1] = 1e7;

    let (sel_i, sel_j) = block_align(&crossscore, params.gap.open);

    let mut anchors: Vec<(usize, usize)> = Vec::new();
    for (&si, &sj) in sel_i.iter().zip(sel_j.iter()) {
        if si == 0 || si == size - 1 { continue; }
        if sj == 0 || sj == size - 1 { continue; }
        let r1 = si - 1;
        let r2 = sj - 1;
        if r1 >= nseg || r2 >= nseg { continue; }
        let seg_idx = sort1[r1];
        if rank2[seg_idx] != r2 { continue; }
        let seg = &all[seg_idx];
        if seg.center1 < n && seg.center2 < m {
            anchors.push((seg.center1, seg.center2));
        }
    }

    if anchors.is_empty() {
        None
    } else {
        Some(anchors)
    }
}

/// Perform FFT-accelerated profile alignment.
///
/// 1. Converts profiles to per-residue-type complex vectors.
/// 2. Uses multi-channel FFT cross-correlation to find the best lags.
/// 3. Detects alignable segments at each candidate lag.
/// 4. Selects optimal non-overlapping anchor pairs via DP (block_align).
/// 5. Runs full DP alignment within each segment.
///
/// Falls back to direct profile DP if no good anchors are found.
pub fn fft_profile_align(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    params: &FftAlignParams,
) -> Alignment {
    if prof1.length == 0 || prof2.length == 0 {
        return Alignment {
            seq1: Vec::new(),
            seq2: Vec::new(),
            score: 0.0,
            operations: Vec::new(),
        };
    }

    match find_fft_anchors(prof1, prof2, matrix, params) {
        Some(anchors) => align_with_anchors(prof1, prof2, matrix, &params.gap, &anchors),
        None => profile_align(prof1, prof2, matrix, &params.gap, params.head_gap, params.tail_gap),
    }
}

/// Convert profile frequencies to per-channel complex vectors for FFT.
fn profile_to_channels(
    prof: &Profile,
    num_channels: usize,
    fft_size: usize,
) -> Vec<Vec<num_complex::Complex64>> {
    let mut channels =
        vec![vec![num_complex::Complex64::new(0.0, 0.0); fft_size]; num_channels];
    for pos in 0..prof.length.min(fft_size) {
        for ch in 0..num_channels.min(prof.nalphabets) {
            channels[ch][pos] = num_complex::Complex64::new(prof.freqs[pos][ch], 0.0);
        }
    }
    channels
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

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_setup() -> (Vec<Vec<i32>>, [u8; 256], usize) {
        let mut mtx = vec![vec![-100i32; 5]; 5];
        for i in 0..4 {
            mtx[i][i] = 100;
        }
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
            num_channels: 4,
        };
        let aln = fft_profile_align(&prof, &prof, &mtx, &params);
        assert!(aln.score > 0.0);
        assert!(!aln.operations.is_empty());
    }

    #[test]
    fn fft_align_falls_back_on_short() {
        let (mtx, map, nalpha) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGT"];
        let prof = Profile::from_aligned(&seqs, &[1.0], &map, nalpha);
        let params = FftAlignParams::protein();
        let aln = fft_profile_align(&prof, &prof, &mtx, &params);
        assert!((aln.score - 400.0).abs() < 1e-6);
    }
}
