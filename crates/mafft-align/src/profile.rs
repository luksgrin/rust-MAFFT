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
            // Frequencies: C's cpmx_calc_new maps ALL characters through amino_n,
            // including '-' (index 24) and '.' (index 23). Gap characters participate
            // in the composition probability matrix and thus in match_score computation.
            for (pos, &ch) in seq.iter().enumerate() {
                if pos >= length { break; }
                let is_gap = ch == b'-' || ch == b'.';
                if is_gap {
                    gap_freq[pos] += w;
                }
                // Map ALL characters (including gaps) into freqs, matching C's cpmx_calc_new
                let idx = amino_map[ch as usize] as usize;
                if idx < nalphabets {
                    freqs[pos][idx] += w;
                }
            }

            // Gap opening count: C's st_OpeningGapCount (mltaln9.c:12837)
            // ogcp[i] counts non-gap→gap transitions at position i.
            // gc starts as 0 (assumes non-gap before position 0).
            {
                let mut gc = false; // gc = 0 in C
                for pos in 0..length {
                    let gb = gc;
                    gc = seq.get(pos).map_or(true, |&c| c == b'-' || c == b'.');
                    if !gb && gc {
                        opening_count[pos] += w;
                    }
                }
            }

            // Gap closing count: C's st_FinalGapCount (mltaln9.c:12880)
            // fgcp[i] counts gap→non-gap transitions where gap is at position i
            // and non-gap is at position i+1. gc starts as seq[0].
            {
                let mut gc = seq.first().map_or(true, |&c| c == b'-' || c == b'.');
                for pos in 0..length {
                    let gb = gc;
                    gc = seq.get(pos + 1).map_or(false, |&c| c == b'-' || c == b'.');
                    // C: gc = 0 at tail (assumes non-gap after last position)
                    if gb && !gc {
                        closing_count[pos] += w;
                    }
                }
            }
        }

        // C: nongap_freq = 1.0 - gap_freq (no clamping).
        // C's convention is that weights sum to 1.0 (normalized before passing),
        // but we preserve C's exact formula even with unnormalized weights
        // to match cpmx_calc_new + gapcountf + st_*GapCount behavior exactly.
        let nongap_freq: Vec<f64> = gap_freq.iter().map(|&g| 1.0 - g).collect();

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
///
/// Segment gap handling matches C's Falign (lines 686-687):
///   - First segment:       headgp = outgap (0), tailgp = 1
///   - Intermediate:        headgp = 1,          tailgp = 1
///   - Last segment:        headgp = 1,          tailgp = outgap (0)
///
/// With outgap=0 (always set by the mafft script via `-O`):
///   first segment gets head_gap=false, last gets tail_gap=false,
///   and all intermediate segments get head_gap=true, tail_gap=true.
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

    // Collect segment boundaries: each anchor defines a break point.
    // Segments are: [0..anchor0), anchor0, [anchor0+1..anchor1), anchor1, ... [lastanchor+1..end)
    // C counts: segment 0 = before first anchor, segment count-2 = after last anchor.
    let num_anchors = anchors.len();

    for (anchor_idx, &(a1, a2)) in anchors.iter().enumerate() {
        if a1 > p1 || a2 > p2 {
            let sub1 = prof1.sub_profile(p1, a1);
            let sub2 = prof2.sub_profile(p2, a2);
            let is_first = anchor_idx == 0;
            // headgp: first segment gets outgap (false), others get 1 (true).
            // tailgp: all pre-anchor segments are followed by an anchor, so tailgp=1 (true).
            let seg_aln = profile_align(&sub1, &sub2, matrix, gap, !is_first, true);
            total_score += seg_aln.score;
            all_ops.extend(seg_aln.operations);
        }
        if a1 < prof1.length && a2 < prof2.length {
            total_score += prof1.match_score(a1, prof2, a2, matrix);
            all_ops.push(AlignOp::Match);
            p1 = a1 + 1;
            p2 = a2 + 1;
        }
    }

    // Trailing segment (after last anchor): headgp=1, tailgp=outgap (false).
    if p1 < prof1.length || p2 < prof2.length {
        let sub1 = prof1.sub_profile(p1, prof1.length);
        let sub2 = prof2.sub_profile(p2, prof2.length);
        let is_only_segment = num_anchors == 0;
        // headgp: if there are anchors, this follows an anchor → headgp=1.
        //         if there are no anchors, this is the first (and only) segment → headgp=outgap (false).
        // tailgp: last segment → tailgp=outgap (false).
        let seg_aln = profile_align(&sub1, &sub2, matrix, gap, !is_only_segment, false);
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
    profile_align_imp(prof1, prof2, matrix, gap, head_gap, tail_gap, None)
}

/// Profile alignment with optional per-cell importance bonuses.
///
/// `impmtx[i][j]` is added to the match score at DP cell `(i, j)`. Mirrors
/// C's `imp_match_out_vead` calls in `Salignmm.c::A__align` (lines 1700-1849)
/// where `currentw[j] += impmtx[i][j]` is applied before each row's cell
/// updates. Used by L-INS-i / E-INS-i to weight DP cells by local-homology
/// importance.
pub fn profile_align_imp(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<i32>],
    gap: &GapModel,
    head_gap: bool,
    tail_gap: bool,
    impmtx: Option<&[Vec<f64>]>,
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
    // Uses C's exact rolling-row scheme, batch match_calc with sparse
    // dot product, and pointer-order-matching float accumulation.

    let hgf1: f64 = 1.0;
    let hgf2: f64 = 1.0;
    let gf1_0 = prof1.nongap_freq.first().copied().unwrap_or(1.0);
    let gf2_0 = prof2.nongap_freq.first().copied().unwrap_or(1.0);
    let nalpha = prof1.nalphabets.min(prof2.nalphabets).min(matrix.len());

    // Build sparse representation of prof2 (C's cpmxpd/cpmxpdn, lines 196-212).
    // For each position j, store only non-zero (alphabet_index, frequency) pairs.
    let cpmx2_sparse: Vec<Vec<(usize, f64)>> = (0..m).map(|j| {
        let mut entries = Vec::new();
        for l in 0..nalpha {
            let v = prof2.freqs[j][l];
            if v != 0.0 {
                entries.push((l, v));
            }
        }
        entries
    }).collect();

    // Batch match_calc: compute scarr once per row, dot with sparse prof2.
    // C's exact loop order (lines 216-223):
    //   for l in 0..nalphabets:
    //     scarr[l] = 0
    //     for j in 0..nalphabets:
    //       scarr[l] += matrix[j][l] * cpmx1[j][position]
    // C's match_calc writes to output[0..lgth2-1] (0-based, no +1 shift).
    // At the tail row (row_pos == n), C's over-allocated arrays give zeros.
    let match_calc_row = |row_pos: usize, output: &mut [f64]| {
        if row_pos >= n {
            // Out-of-bounds: C's calloc'd arrays produce zeros here
            for j in 0..m {
                output[j] = 0.0;
            }
            return;
        }
        let mut scarr = vec![0.0f64; nalpha];
        for l in 0..nalpha {
            scarr[l] = 0.0;
            for j in 0..nalpha {
                scarr[l] += matrix[j][l] as f64 * prof1.freqs[row_pos][j];
            }
        }
        for j in 0..m {
            output[j] = 0.0;
            for &(k, v) in &cpmx2_sparse[j] {
                output[j] += scarr[k] * v;
            }
        }
    };

    // h matrix (need full for traceback) and ijp
    let mut h = vec![vec![0.0f64; m + 1]; n + 1];
    let mut ijp = vec![vec![0i32; m + 1]; n + 1];

    // initverticalw (C line 776): match_calc with prof2 pos 0 vs all prof1 positions.
    // C calls match_calc with swapped profiles: cpmx2pt first, cpmx1pt second.
    let cpmx1_sparse: Vec<Vec<(usize, f64)>> = (0..n).map(|i| {
        let mut entries = Vec::new();
        for l in 0..nalpha {
            let v = prof1.freqs[i][l];
            if v != 0.0 {
                entries.push((l, v));
            }
        }
        entries
    }).collect();

    // initverticalw (C line 776): match_calc(cpmx2pt, cpmx1pt, 0, lgth1, initverticalw)
    // C fills initverticalw[0..lgth1-1] (0-based), then adds gap to [1..lgth1].
    // Result: [0] = pure match score (no gap), [1..n-1] = match + gap, [n] = 0 + gap.
    let mut initverticalw = vec![0.0f64; n + 1];
    {
        let mut scarr = vec![0.0f64; nalpha];
        for l in 0..nalpha {
            scarr[l] = 0.0;
            for j in 0..nalpha {
                scarr[l] += matrix[j][l] as f64 * prof2.freqs[0][j];
            }
        }
        for i in 0..n {
            initverticalw[i] = 0.0;
            for &(k, v) in &cpmx1_sparse[i] {
                initverticalw[i] += scarr[k] * v;
            }
        }
    }
    // C: imp_match_out_vead_tate(initverticalw, 0, lgth1) — add impmtx[i][0] to initverticalw[i].
    if let Some(imp) = impmtx {
        for i in 0..n {
            initverticalw[i] += imp[i][0];
        }
    }
    if head_gap {
        for i in 1..=n {
            initverticalw[i] += ogcp1[0] * hgf2 + fgcp1[i - 1] * gf2_0;
        }
    }

    // currentw (C line 779): match_calc with prof1 pos 0 vs all prof2 positions.
    // C fills currentw[0..lgth2-1] (0-based), then adds gap to [1..lgth2].
    // Result: [0] = pure match score (no gap), [1..m-1] = match + gap, [m] = 0 + gap.
    let mut currentw = vec![0.0f64; m + 1];
    {
        let mut scarr = vec![0.0f64; nalpha];
        for l in 0..nalpha {
            scarr[l] = 0.0;
            for j in 0..nalpha {
                scarr[l] += matrix[j][l] as f64 * prof1.freqs[0][j];
            }
        }
        for j in 0..m {
            currentw[j] = 0.0;
            for &(k, v) in &cpmx2_sparse[j] {
                currentw[j] += scarr[k] * v;
            }
        }
    }
    // C: imp_match_out_vead(currentw, 0, lgth2) — add impmtx[0][j] to currentw[j].
    if let Some(imp) = impmtx {
        for j in 0..m {
            currentw[j] += imp[0][j];
        }
    }
    if head_gap {
        for j in 1..=m {
            currentw[j] += ogcp2[0] * hgf1 + fgcp2[j - 1] * gf1_0;
        }
    }

    for j in 0..=m { h[0][j] = currentw[j]; }
    for i in 0..=n { h[i][0] = initverticalw[i]; }

    let mut mj = vec![f64::NEG_INFINITY; m + 1];
    let mut mpj = vec![0usize; m + 1];
    for j in 1..=m {
        let gf2_jm1 = prof2.nongap_freq.get(j - 1).copied().unwrap_or(0.0);
        mj[j] = currentw[j - 1] + ogcp1[1] * gf2_jm1;
        mpj[j] = 0;
    }

    for i in 0..=n { ijp[i][0] = i as i32 + 1; }
    for j in 0..=m { ijp[0][j] = -(j as i32 + 1); }

    let mut previousw = vec![0.0f64; m + 1];
    let mut lastverticalw = vec![0.0f64; n + 1];
    lastverticalw[0] = currentw[m - 1];

    // C pads gapfreq1pt[lgth1] = 1.0 and gapfreq2pt[lgth2] = 1.0 (tditeration.c's
    // `for(i=0;i<lgth+1;i++) gapfreq[i] = 1.0 - gapfreq[i];` with calloc'd 0 → 1).
    // Our DP needs these boundary values when i==n or j==m.
    let lasti = if tail_gap { n + 1 } else { n };
    for i in 1..lasti {
        std::mem::swap(&mut previousw, &mut currentw);
        previousw[0] = initverticalw[i - 1];

        // Batch match_calc for row i (C line 816: ist+i)
        // Fills currentw[0..m-1]. Position m doesn't exist in prof2.
        match_calc_row(i, &mut currentw);
        // C: imp_match_out_vead(currentw, i, lgth2) — add impmtx[i][:] (Salignmm.c:1848).
        if let Some(imp) = impmtx {
            if i < imp.len() {
                let row = &imp[i];
                for j in 0..m.min(row.len()) {
                    currentw[j] += row[j];
                }
            }
        }
        currentw[m] = 0.0; // padding: no substitution score beyond prof2
        currentw[0] = initverticalw[i];

        let gf1_im1 = prof1.nongap_freq[i - 1]; // i-1 in 0..n-1, always valid
        let mut mi = previousw[0] + ogcp2[1] * gf1_im1;
        let mut mpi: usize = 0;

        for j in 1..=m {
            // Out-of-bounds positions in nongap_freq correspond to C's padded
            // gapfreq1pt[lgth1] / gapfreq2pt[lgth2] = 1.0.
            let gf1_i = if i < n { prof1.nongap_freq[i] } else { 1.0 };
            let gf1_im1 = prof1.nongap_freq[i - 1];
            let gf2_j = if j < m { prof2.nongap_freq[j] } else { 1.0 };
            let gf2_jm1 = prof2.nongap_freq[j - 1];

            let mut wm = previousw[j - 1];
            ijp[i][j] = 0;

            let g = mi + fgcp2[j - 1] * gf1_i;
            if g > wm {
                wm = g;
                ijp[i][j] = -(j as i32 - mpi as i32);
            }

            let g = previousw[j - 1] + ogcp2[j] * gf1_im1;
            if g >= mi {
                mi = g;
                mpi = j - 1;
            }

            let g = mj[j] + fgcp1[i - 1] * gf2_j;
            if g > wm {
                wm = g;
                ijp[i][j] = i as i32 - mpj[j] as i32;
            }

            let g = previousw[j - 1] + ogcp1[i] * gf2_jm1;
            if g >= mj[j] {
                mj[j] = g;
                mpj[j] = i - 1;
            }

            currentw[j] += wm;
            h[i][j] = currentw[j];
        }
        lastverticalw[i] = currentw[m - 1];
    }

    // Tail gap handling  (C lines 512-536)
    // For tail_gap=false, the main DP loop stops at row n-1 (lasti = n).
    // Find the best ending position in the last computed row or last column.
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
        // Set h[n][m] to the best score found (for the score field).
        h[n][m] = wm;
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

/// Pairwise alignment for single-sequence pairs.
///
/// Ports C's `G__align11()` from Galign11.c. Uses character-indexed scoring
/// matrix and flat gap penalties (no position-specific ogcp/fgcp weighting).
/// C uses this instead of MSalignmm when both groups have exactly 1 sequence.
pub fn pairwise_align11(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    penalty: f64,
    head_gap: bool,
    tail_gap: bool,
) -> Alignment {
    let n = seq1.len();
    let m = seq2.len();
    if n == 0 || m == 0 {
        return Alignment { seq1: Vec::new(), seq2: Vec::new(), score: 0.0, operations: Vec::new() };
    }

    let nalpha = matrix.len();

    // Build amino_dynamicmtx[char][char] (C line 1128-1129)
    // For direct character lookup without going through amino_map during DP
    let mut amino_mtx = vec![vec![0.0f64; 256]; 256];
    // We need the reverse map: index -> char. Build from amino_map.
    // Actually, just use amino_map inline during match_calc.

    // match_calc_mtx: score = amino_dynamicmtx[seq1[i]][seq2[j]]
    // Which equals matrix[amino_map[seq1[i]]][amino_map[seq2[j]]]
    let score_pair = |c1: u8, c2: u8| -> f64 {
        let i = amino_map[c1 as usize] as usize;
        let j = amino_map[c2 as usize] as usize;
        if i < nalpha && j < nalpha { matrix[i][j] as f64 } else { 0.0 }
    };

    // initverticalw: match_calc_mtx(seq2, seq1, 0, lgth1)
    // C fills initverticalw[0..lgth1-1] (0-based), then adds gap to [1..lgth1].
    // Result: [0] = pure match score (no gap), [1..n-1] = match + gap, [n] = 0 + gap.
    let mut initverticalw = vec![0.0f64; n + 1];
    for i in 0..n {
        initverticalw[i] = score_pair(seq1[i], seq2[0]);
    }
    if head_gap {
        for i in 1..=n {
            initverticalw[i] += penalty; // flat penalty, C line 1180
        }
    }

    // currentw: match_calc_mtx(seq1, seq2, 0, lgth2)
    // C fills currentw[0..lgth2-1] (0-based), then adds gap to [1..lgth2].
    let mut currentw = vec![0.0f64; m + 1];
    for j in 0..m {
        currentw[j] = score_pair(seq1[0], seq2[j]);
    }
    if head_gap {
        for j in 1..=m {
            currentw[j] += penalty; // flat penalty, C line 1188
        }
    }

    let mut h = vec![vec![0.0f64; m + 1]; n + 1];
    let mut ijp = vec![vec![0i32; m + 1]; n + 1];

    for j in 0..=m { h[0][j] = currentw[j]; }
    for i in 0..=n { h[i][0] = initverticalw[i]; }

    // m[j] = currentw[j-1], NO ogcp multiplication (C line 1218)
    let mut mj = vec![f64::NEG_INFINITY; m + 1];
    let mut mpj = vec![0usize; m + 1];
    for j in 1..=m {
        mj[j] = currentw[j - 1];
        mpj[j] = 0;
    }

    for i in 0..=n { ijp[i][0] = i as i32 + 1; }
    for j in 0..=m { ijp[0][j] = -(j as i32 + 1); }

    let mut previousw = vec![0.0f64; m + 1];
    let mut lastverticalw = vec![0.0f64; n + 1];
    if m > 0 { lastverticalw[0] = currentw[m - 1]; }

    let lasti = if tail_gap { n + 1 } else { n };
    for i in 1..lasti {
        std::mem::swap(&mut previousw, &mut currentw);
        previousw[0] = initverticalw[i - 1];

        // match_calc_mtx for row i (C line 1252: uses index i, not i-1)
        // C fills currentw[0..lgth2-1]. Position m is padding (0).
        if i < n {
            for j in 0..m {
                currentw[j] = score_pair(seq1[i], seq2[j]);
            }
        } else {
            // Out-of-bounds tail row: C's over-allocated arrays give zeros
            for j in 0..m {
                currentw[j] = 0.0;
            }
        }
        currentw[m] = 0.0; // padding: no substitution score beyond seq2
        currentw[0] = initverticalw[i];

        // mi = previousw[0], NO ogcp2 multiplication (C line 1275)
        let mut mi = previousw[0];
        let mut mpi: usize = 0;

        for j in 1..=m {
            let mut wm = previousw[j - 1];
            ijp[i][j] = 0;

            // Deletion: mi + fpenalty (flat, C line 1306)
            let g = mi + penalty;
            if g > wm {
                wm = g;
                ijp[i][j] = -(j as i32 - mpi as i32);
            }
            // Open from diagonal (C line 1311)
            let g = previousw[j - 1];
            if g >= mi {
                mi = g;
                mpi = j - 1;
            }

            // Insertion: m[j] + fpenalty (flat, C line 1325)
            let g = mj[j] + penalty;
            if g > wm {
                wm = g;
                ijp[i][j] = i as i32 - mpj[j] as i32;
            }
            // Open from diagonal (C line 1330)
            let g = previousw[j - 1];
            if g >= mj[j] {
                mj[j] = g;
                mpj[j] = i - 1;
            }

            currentw[j] += wm;
            h[i][j] = currentw[j];
        }
        if m > 0 { lastverticalw[i] = currentw[m - 1]; }
    }

    // Tail gap handling
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
        // Set h[n][m] to the best score found (for the score field).
        h[n][m] = wm;
    }

    // Traceback (same ijp format as MSalignmm)
    let mut gaptable1 = Vec::new();
    let mut gaptable2 = Vec::new();
    let (mut iin, mut jin) = (n as i32, m as i32);
    let klim = n + m;
    let mut k = 0;
    while k <= klim {
        let v = ijp[iin as usize][jin as usize];
        let (ifi, jfi): (i32, i32) = if v < 0 {
            (iin - 1, jin + v)
        } else if v > 0 {
            (iin - v, jin - 1)
        } else {
            (iin - 1, jin - 1)
        };
        let mut l = iin - ifi;
        while l > 1 { gaptable1.push(b'o'); gaptable2.push(b'-'); k += 1; l -= 1; }
        let mut l = jin - jfi;
        while l > 1 { gaptable1.push(b'-'); gaptable2.push(b'o'); k += 1; l -= 1; }
        if iin <= 0 || jin <= 0 { break; }
        gaptable1.push(b'o'); gaptable2.push(b'o'); k += 1;
        iin = ifi; jin = jfi;
    }

    gaptable1.reverse();
    gaptable2.reverse();
    let mut ops = Vec::with_capacity(gaptable1.len());
    for k in 0..gaptable1.len() {
        match (gaptable1[k], gaptable2[k]) {
            (b'o', b'o') => ops.push(AlignOp::Match),
            (b'o', b'-') => ops.push(AlignOp::Delete),
            (b'-', b'o') => ops.push(AlignOp::Insert),
            _ => {}
        }
    }

    Alignment { seq1: Vec::new(), seq2: Vec::new(), score: h[n][m], operations: ops }
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
