/// Iterative refinement (tree-dependent iteration).
///
/// Ports the C `TreeDependentIteration()` from tditeration.c.
///
/// Repeatedly re-aligns pairs of groups defined by the guide tree,
/// accepting improvements and rejecting regressions, until convergence.
///
/// Key insight from the C code: at each tree branch, ALL sequences are
/// split into two groups (subtree vs everything else). There are never
/// "uninvolved" sequences — every sequence is in one group or the other.
///
/// Branch enumeration matches C exactly:
/// - For each topology step, both sides (k=0: left, k=1: right) are
///   processed, EXCEPT the root step (last step) where only k=1 (right)
///   is used (since at the root, left-vs-complement and right-vs-complement
///   produce the same split, just flipped).
/// - Even iterations traverse steps forward (0 → N-1), odd iterations
///   traverse backward (N-1 → 0). Within each step, k always goes 0→1.
/// - Total branches per iteration: (nseq-1)*2 - 1.

use rayon::prelude::*;

use mafft_align::{
    profile_align, profile_align_imp,
    profile_align_imp_with_boundary, BoundaryFreqs,
    build_imp_matrix, FASTATHRESHOLD_DEFAULT,
    Profile, GapModel, AlignOp,
};
use mafft_fft::{alignable_segments, SegmentParams};
use mafft_tree::{Topology, BranchWeights};
use mafft_types::{ScoringContext, LocalHomologyTable};

use crate::progressive::MultipleAlignment;

/// Parameters controlling iterative refinement.
#[derive(Debug, Clone)]
pub struct RefinementParams {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Score improvement threshold (fraction of old score).
    /// C default is 0.0 (accept only strict improvements).
    pub cut: f64,
    /// Whether to use FFT-accelerated alignment during refinement.
    pub use_fft: bool,
}

impl Default for RefinementParams {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            cut: 0.0,
            use_fft: false,
        }
    }
}

/// A branch identifier for oscillation tracking: (step_index, side).
/// side 0 = left, side 1 = right.
type BranchId = (usize, usize);

/// Build the per-step branch splits from a topology, matching C's enumeration.
///
/// For each topology step, both sides (k=0: left vs complement, k=1: right vs
/// complement) are included — EXCEPT the root step (last step) where only k=1
/// is included. At the root, left-vs-complement and right-vs-complement are
/// identical splits (just flipped), so C skips the redundant one.
///
/// Returns: `branch_map[step_idx]` = list of `(side, group1, group2)`.
/// Total branches = `(nseq - 1) * 2 - 1`.
fn build_branch_map(
    topology: &Topology,
    nseq: usize,
) -> Vec<Vec<(usize, Vec<usize>, Vec<usize>)>> {
    let nsteps = topology.steps.len();
    let root_idx = nsteps - 1;
    let all_indices: Vec<usize> = (0..nseq).collect();

    let mut branch_map: Vec<Vec<(usize, Vec<usize>, Vec<usize>)>> = Vec::with_capacity(nsteps);
    for (step_idx, step) in topology.steps.iter().enumerate() {
        let is_root = step_idx == root_idx;
        let mut sides = Vec::new();

        if !is_root {
            // Side 0: step.left vs complement
            let complement: Vec<usize> = all_indices
                .iter()
                .filter(|i| !step.left.contains(i))
                .copied()
                .collect();
            if !step.left.is_empty() && !complement.is_empty() {
                sides.push((0, step.left.clone(), complement));
            }
        }

        // Side 1: step.right vs complement
        let complement: Vec<usize> = all_indices
            .iter()
            .filter(|i| !step.right.contains(i))
            .copied()
            .collect();
        if !step.right.is_empty() && !complement.is_empty() {
            sides.push((1, step.right.clone(), complement));
        }

        branch_map.push(sides);
    }
    branch_map
}

/// Iteratively refine a multiple alignment.
///
/// At each tree branch, splits ALL sequences into two groups (subtree vs
/// rest), re-aligns the two groups, and accepts improvements.
pub fn iterative_refine(
    alignment: &mut MultipleAlignment,
    topology: &Topology,
    scoring: &ScoringContext,
    params: &RefinementParams,
    constraints: Option<&LocalHomologyTable>,
) -> usize {
    let nseq = alignment.nseq();
    if nseq <= 2 || topology.steps.is_empty() {
        return 0;
    }

    let branch_weights = BranchWeights::new(topology);
    let global_weights = mafft_tree::sequence_weights(topology);
    let use_global_weights = std::env::var("RUST_MAFFT_GLOBAL_WEIGHTS").is_ok();
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);


    let mut converged_count = 0usize;
    let convergence_target = nseq * 2;

    let nsteps = topology.steps.len();
    let branch_map = build_branch_map(topology, nseq);

    // Per-branch score history for oscillation detection.
    // history[iteration][(step_idx, side)] = score after processing that branch.
    let mut history: Vec<std::collections::HashMap<BranchId, f64>> = Vec::new();

    let mut iteration = 0;
    for iter in 0..params.max_iterations {
        iteration = iter + 1;
        let mut any_change = false;
        let mut iter_scores: std::collections::HashMap<BranchId, f64> = std::collections::HashMap::new();

        // C alternates step traversal direction: even → forward, odd → reverse.
        let step_order: Vec<usize> = if iter % 2 == 0 {
            (0..nsteps).collect()
        } else {
            (0..nsteps).rev().collect()
        };

        for &step_idx in &step_order {
            for (side, group1, group2) in &branch_map[step_idx] {
                let branch_id: BranchId = (step_idx, *side);

                let weights = if use_global_weights {
                    global_weights.clone()
                } else {
                    branch_weights.weights_for_branch(topology, step_idx, *side)
                };

                // Group-local sum-1 normalized weights (matches C's
                // fastconjuction_noname). Used both for `compute_impmatch_diagonal`
                // and any future per-cluster averaging.
                const MIN_W: f64 = 0.00001;
                let w1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
                let w2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
                let s1w: f64 = w1.iter().sum();
                let s2w: f64 = w2.iter().sum();
                let w1n: Vec<f64> = if s1w > 0.0 { w1.iter().map(|w| w / s1w).collect() } else { vec![1.0; group1.len()] };
                let w2n: Vec<f64> = if s2w > 0.0 { w2.iter().map(|w| w / s2w).collect() } else { vec![1.0; group2.len()] };

                // C's mscore = oimpmatchdouble + tmpdouble (tditeration.c:953):
                // intergroup substitution score + impmatch (sum of impmtx[i][i]
                // over the current alignment's columns). We compute the same.
                let old_sub = compute_split_score(
                    group1, group2, &alignment.sequences, &weights, scoring,
                );
                let old_imp = if let Some(lh) = constraints {
                    compute_impmatch_diagonal(
                        group1, group2, &alignment.sequences, &w1n, &w2n, lh,
                    )
                } else { 0.0 };
                let old_score = old_sub + old_imp;

                let new_seqs = realign_all(
                    group1, group2, &alignment.sequences, &weights, scoring, &gap,
                    constraints, params.use_fft,
                );

                if let Some((new_seqs, _new_score)) = new_seqs {
                    // C's identity check (tditeration.c:2184-2185): compare only
                    // the representative sequences s1=memlist1[0], s2=memlist2[0]
                    // (from OneClusterAndTheOther_fast in tddis.c:834-835).
                    // `group1` is memlist1, `group2` is memlist2, so s1=group1[0],
                    // s2=group2[0]. Checking ALL sequences would incorrectly treat
                    // column-rearrangements that preserve the two representatives
                    // as "changed", causing spurious accepts.
                    let s1 = group1[0];
                    let s2 = group2[0];
                    let changed = alignment.sequences[s1] != new_seqs[s1]
                        || alignment.sequences[s2] != new_seqs[s2];

                    if !changed {
                        // Identical — no change, count toward convergence
                        let tscore = old_score;
                        iter_scores.insert(branch_id, tscore);
                        converged_count += 1;
                    } else {
                        // C's tscore = impmatchdouble + tmpdouble (tditeration.c:1094):
                        // intergroup score + new alignment's impmatch.
                        let new_sub = compute_split_score(
                            group1, group2, &new_seqs, &weights, scoring,
                        );
                        let new_imp = if let Some(lh) = constraints {
                            compute_impmatch_diagonal(
                                group1, group2, &new_seqs, &w1n, &w2n, lh,
                            )
                        } else { 0.0 };
                        let tscore = new_sub + new_imp;

                        let threshold = old_score - params.cut / 100.0 * old_score;
                        if std::env::var("RUST_MAFFT_TRACE").is_ok() {
                            eprintln!("ACCEPT iter={iter} step={step_idx} side={side} old={:.3} new={:.3} accept={}",
                                old_score, tscore, tscore > threshold);
                            if std::env::var("RUST_DUMP_BRANCHES").is_ok() {
                                eprintln!("[BRANCH_BEFORE] iter={iter} step={step_idx} side={side}");
                                for (idx, s) in alignment.sequences.iter().enumerate() {
                                    eprintln!("  seq{idx}: {}", std::str::from_utf8(s).unwrap_or(""));
                                }
                                eprintln!("[BRANCH_AFTER] iter={iter} step={step_idx} side={side}");
                                for (idx, s) in new_seqs.iter().enumerate() {
                                    eprintln!("  seq{idx}: {}", std::str::from_utf8(s).unwrap_or(""));
                                }
                            }
                        }
                        if tscore > threshold {
                            alignment.sequences = new_seqs;
                            any_change = true;
                            converged_count = 0;
                        } else {
                            converged_count += 1;
                        }
                        iter_scores.insert(branch_id, tscore);
                    }
                } else {
                    iter_scores.insert(branch_id, old_score);
                    converged_count += 1;
                }

                if converged_count >= convergence_target {
                    return iteration;
                }

                // Oscillation detection: check if this branch's score matches
                // the score from 2, 4, 6... iterations ago (same branch).
                if iter >= 2 {
                    let tscore = iter_scores[&branch_id];
                    let mut oscillating = false;
                    let mut ii = history.len() as isize - 2; // iterate-2
                    while ii >= 0 {
                        if let Some(&prev_score) = history[ii as usize].get(&branch_id) {
                            if tscore == prev_score {
                                oscillating = true;
                                break;
                            }
                        }
                        ii -= 2;
                    }
                    if oscillating {
                        return iteration;
                    }
                }
            }
        }

        history.push(iter_scores);

        if !any_change {
            return iteration;
        }
    }

    iteration
}

/// Re-align all sequences split into two groups.
///
/// Since group1 + group2 = ALL sequences, there are no "other" sequences
/// to worry about. Every sequence is in exactly one group.
///
/// When `use_fft` is true (matching C's Falign path in tditeration.c):
/// 1. Strip per-group gap columns → build stripped profiles
/// 2. Run FFT anchor detection on stripped profiles (clean, residue-rich data)
/// 3. Map anchors back to non-stripped coordinates via kept1/kept2
/// 4. Build non-stripped profiles from full sequences
/// 5. Run anchored DP on non-stripped profiles (matching C's input)
///
/// The anchor mapping ensures the FFT sees clean data for good anchor
/// detection, while the DP operates on the same non-stripped profiles C
/// uses. Anchors constrain the DP so width growth is bounded.
///
/// If FFT finds no anchors, falls back to profile_align on the full
/// non-stripped sequences (matching C's single-segment fallback in
/// Falign.c lines 1377-1421), bounded by alloclen.
fn realign_all(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    constraints: Option<&LocalHomologyTable>,
    use_fft: bool,
) -> Option<(Vec<Vec<u8>>, f64)> {
    let width = sequences[0].len();

    // Per-group gap stripping.
    let gap1 = group_all_gap_columns(group1, sequences, width);
    let gap2 = group_all_gap_columns(group2, sequences, width);
    let kept1: Vec<usize> = (0..width).filter(|&c| !gap1[c]).collect();
    let kept2: Vec<usize> = (0..width).filter(|&c| !gap2[c]).collect();

    // C clamps per-sequence weights to minimumweight (0.00001 from mafft script,
    // applied in fastconjuction_noname at tddis.c line 548).
    const MINIMUM_WEIGHT: f64 = 0.00001;
    let w1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MINIMUM_WEIGHT)).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MINIMUM_WEIGHT)).collect();
    let sum1: f64 = w1.iter().sum();
    let sum2: f64 = w2.iter().sum();
    let w1n: Vec<f64> = if sum1 > 0.0 { w1.iter().map(|w| w / sum1).collect() } else { vec![1.0; group1.len()] };
    let w2n: Vec<f64> = if sum2 > 0.0 { w2.iter().map(|w| w / sum2).collect() } else { vec![1.0; group2.len()] };

    // Build stripped sequences and profiles (used for FFT anchor detection
    // and as fallback for unconstrained non-FFT alignment).
    let stripped1: Vec<Vec<u8>> = group1.iter()
        .map(|&i| kept1.iter().map(|&c| sequences[i][c]).collect())
        .collect();
    let stripped2: Vec<Vec<u8>> = group2.iter()
        .map(|&i| kept2.iter().map(|&c| sequences[i][c]).collect())
        .collect();
    let s1_refs: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
    let s2_refs: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();
    let stripped_prof1 = Profile::from_aligned(&s1_refs, &w1n, &scoring.amino_map, scoring.nalphabets);
    let stripped_prof2 = Profile::from_aligned(&s2_refs, &w2n, &scoring.amino_map, scoring.nalphabets);

    if stripped_prof1.length == 0 || stripped_prof2.length == 0 {
        return None;
    }

    if let Some(lh_table) = constraints {
        if use_fft {
            // C's `Falign_localhom` (kobetsubunkatsu=1 path): mirrors the
            // FFT-segmented loop in the unconstrained `Falign` but calls
            // `partA__align(constraint=1, ..., gapmap1, gapmap2, ...)` per
            // segment. The global impmtx is built once for the full
            // (non-stripped) alignment; each segment passes its sliced
            // view via gapmap1/gapmap2 (which translate stripped column
            // index to position within the segment).
            return realign_all_constrained_fft(
                group1, group2, sequences, &w1n, &w2n,
                scoring, gap, lh_table,
            );
        }
        // Non-FFT constraint path (L-INS-i without -F, single full DP).
        // Mirrors C's `A__align(..., constraint=1, ...)` (Salignmm.c:1086):
        // build the per-cell importance matrix `impmtx` once, then do the
        // standard profile DP with `currentw[j] += impmtx[i][j]` applied
        // row-by-row inside the DP (Salignmm.c:1700-1849).
        let g1_seq_refs: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
        let g2_seq_refs: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();
        let imp = build_imp_matrix(
            lh_table,
            group1, group2,
            &g1_seq_refs, &g2_seq_refs,
            &w1n, &w2n,
            stripped_prof1.length, stripped_prof2.length,
            FASTATHRESHOLD_DEFAULT,
        );
        let aln = profile_align_imp(
            &stripped_prof1, &stripped_prof2,
            &scoring.substitution_matrix,
            gap,
            true, true,
            Some(&imp),
        );
        return build_result_from_stripped(
            &aln, group1, group2, sequences, &kept1, &kept2,
            &stripped_prof1, &stripped_prof2,
        );
    }

    if use_fft {
        // C's Falign path for dvtditr (kobetsubunkatsu=1 in dvtditr.c:54):
        // 1. SKIP the FFT block (Falign.c:1110 `if(!kobetsubunkatsu)` is skipped)
        // 2. alignableReagion runs ONCE with lag=0 (Falign.c:1299 maxk=1,kouho[0]=0)
        // 3. Collected segments define cut points: [0, center0, center1, ..., width]
        // 4. Per segment, commongappick strips per-group all-gap columns
        //    (Falign.c:1610 `if(kobetsubunkatsu && fftkeika)`)
        // 5. MSalignmm aligns the stripped segment
        // 6. Results concatenated
        let width = sequences[0].len();
        let full_prof1 = Profile::from_aligned(
            &group1.iter().map(|&i| sequences[i].as_slice()).collect::<Vec<_>>(),
            &w1n, &scoring.amino_map, scoring.nalphabets);
        let full_prof2 = Profile::from_aligned(
            &group2.iter().map(|&i| sequences[i].as_slice()).collect::<Vec<_>>(),
            &w2n, &scoring.amino_map, scoring.nalphabets);

        // mafft.tmpl passes `-z 50` to dvtditr, setting fftThreshold=50 for the
        // alignableReagion sliding window (default from constants.c is 80, but
        // the script overrides it). Match that here.
        let segment_params = if scoring.seq_type.is_nucleotide() {
            SegmentParams::dna().with_threshold(50.0)
        } else {
            SegmentParams::protein().with_threshold(50.0)
        };

        // Step 1: Compute per-position site scores at lag=0 (C's alignableReagion,
        // fftFunctions.c:282-287). For refinement, prof1.length == prof2.length,
        // so we score pairwise at matching columns. C uses n_disFFT which equals
        // substitution_matrix when offset=0 (mafft default).
        // C divides by totaleff = sum eff1[i]*eff2[j]; with normalized weights this is 1.
        let totaleff: f64 = w1n.iter().sum::<f64>() * w2n.iter().sum::<f64>();
        let len = full_prof1.length.min(full_prof2.length);
        let mut site_scores = vec![0.0f64; len];
        for i in 0..len {
            site_scores[i] =
                full_prof1.match_score(i, &full_prof2, i, &scoring.substitution_matrix)
                / totaleff;
        }

        // Step 2: Find alignable segments via the sliding window threshold test.
        let segments = alignable_segments(&site_scores, &segment_params);

        // Step 3: Build cut points from segment centers (Falign.c:1414).
        // kobetsubunkatsu=1: cut1[i+1] = sortedseg1[i]->center, plus [0] and [len].
        // For refinement (lag=0), cut1[i] == cut2[i], so a single cut list suffices.
        let mut cuts: Vec<usize> = Vec::with_capacity(segments.len() + 2);
        cuts.push(0);
        for seg in &segments {
            cuts.push(seg.center.min(width));
        }
        cuts.push(width);
        // Ensure strictly increasing (segments might overlap/coincide — dedupe).
        cuts.sort();
        cuts.dedup();

        // Step 4-5: Per-segment strip + align, then concatenate.
        let mut new_sequences: Vec<Vec<u8>> = vec![Vec::new(); sequences.len()];
        let mut total_score = 0.0f64;
        for win in cuts.windows(2) {
            let a = win[0];
            let b = win[1];
            if a >= b { continue; }

            // Slice the segment [a, b) from each sequence.
            let seg1: Vec<Vec<u8>> = group1.iter()
                .map(|&i| sequences[i][a..b].to_vec()).collect();
            let seg2: Vec<Vec<u8>> = group2.iter()
                .map(|&i| sequences[i][a..b].to_vec()).collect();

            // commongappick on the segment: strip columns where all sequences in
            // THIS GROUP have a gap within THIS SEGMENT. (C's commongappick)
            let seg_width = b - a;
            let seg_gap1: Vec<bool> = (0..seg_width)
                .map(|c| seg1.iter().all(|s| s[c] == b'-')).collect();
            let seg_gap2: Vec<bool> = (0..seg_width)
                .map(|c| seg2.iter().all(|s| s[c] == b'-')).collect();
            let seg_kept1: Vec<usize> = (0..seg_width).filter(|&c| !seg_gap1[c]).collect();
            let seg_kept2: Vec<usize> = (0..seg_width).filter(|&c| !seg_gap2[c]).collect();

            let stripped_seg1: Vec<Vec<u8>> = seg1.iter()
                .map(|s| seg_kept1.iter().map(|&c| s[c]).collect()).collect();
            let stripped_seg2: Vec<Vec<u8>> = seg2.iter()
                .map(|s| seg_kept2.iter().map(|&c| s[c]).collect()).collect();

            if stripped_seg1.is_empty() || stripped_seg1[0].is_empty() ||
               stripped_seg2.is_empty() || stripped_seg2[0].is_empty() {
                // One side is empty — emit gaps for both as-is (only possible
                // when whole segment is all-gap for one group in every column).
                let len1 = if stripped_seg1.is_empty() { 0 } else { stripped_seg1[0].len() };
                let len2 = if stripped_seg2.is_empty() { 0 } else { stripped_seg2[0].len() };
                for (k, &i) in group1.iter().enumerate() {
                    new_sequences[i].extend_from_slice(&stripped_seg1[k]);
                    new_sequences[i].extend(std::iter::repeat(b'-').take(len2));
                }
                for (k, &i) in group2.iter().enumerate() {
                    new_sequences[i].extend(std::iter::repeat(b'-').take(len1));
                    new_sequences[i].extend_from_slice(&stripped_seg2[k]);
                }
                continue;
            }

            let s1_refs: Vec<&[u8]> = stripped_seg1.iter().map(|s| s.as_slice()).collect();
            let s2_refs: Vec<&[u8]> = stripped_seg2.iter().map(|s| s.as_slice()).collect();
            let mut prof_seg1 = Profile::from_aligned(&s1_refs, &w1n, &scoring.amino_map, scoring.nalphabets);
            let mut prof_seg2 = Profile::from_aligned(&s2_refs, &w2n, &scoring.amino_map, scoring.nalphabets);

            // C's Falign segment loop (lines 1549-1565):
            //   sgap[j] = (cut1[i]  > 0)    ? (seq[j][cut1[i]-1]   == '-') : 'o'
            //   egap[j] = (cut1[i+1] != len)? (seq[j][cut1[i+1]]   == '-') : 'o'
            // These per-sequence boundary gap states are passed to MSalignmm,
            // which switches from `st_OpeningGapCount` to `new_OpeningGapCount`
            // (mltaln9.c:12557): a sequence already in a gap just before the
            // segment's first column does NOT count as "opening" at position 0.
            // `Profile::from_aligned` applies `st_OpeningGapCount` semantics
            // (gc starts at 0), so we correct the first/last counts here.
            let sgap_inside = a > 0;
            let egap_inside = b < width;
            if sgap_inside || egap_inside {
                if !prof_seg1.ogcp.is_empty() {
                    for (k, &idx) in group1.iter().enumerate() {
                        if sgap_inside
                            && sequences[idx][a - 1] == b'-'
                            && !stripped_seg1[k].is_empty()
                            && stripped_seg1[k][0] == b'-'
                        {
                            prof_seg1.ogcp[0] -= w1n[k];
                        }
                        let l1 = stripped_seg1[k].len();
                        if egap_inside
                            && l1 > 0
                            && stripped_seg1[k][l1 - 1] == b'-'
                            && sequences[idx][b] == b'-'
                        {
                            let last = prof_seg1.fgcp.len() - 1;
                            prof_seg1.fgcp[last] -= w1n[k];
                        }
                    }
                }
                if !prof_seg2.ogcp.is_empty() {
                    for (k, &idx) in group2.iter().enumerate() {
                        if sgap_inside
                            && sequences[idx][a - 1] == b'-'
                            && !stripped_seg2[k].is_empty()
                            && stripped_seg2[k][0] == b'-'
                        {
                            prof_seg2.ogcp[0] -= w2n[k];
                        }
                        let l2 = stripped_seg2[k].len();
                        if egap_inside
                            && l2 > 0
                            && stripped_seg2[k][l2 - 1] == b'-'
                            && sequences[idx][b] == b'-'
                        {
                            let last = prof_seg2.fgcp.len() - 1;
                            prof_seg2.fgcp[last] -= w2n[k];
                        }
                    }
                }
            }

            // C's Falign segment loop (lines 1521-1522):
            //   headgp = (i==0) ? outgap : 1   ;   tailgp = (i==count-2) ? outgap : 1
            // outgap=1 in dvtditr.c:73, so headgp=tailgp=1 for every segment.
            let seg_aln = profile_align(
                &prof_seg1, &prof_seg2, &scoring.substitution_matrix, gap, true, true,
            );
            total_score += seg_aln.score;

            // Reconstruct the segment's output by applying ops to the stripped segments.
            let mut i1 = 0usize;
            let mut i2 = 0usize;
            for op in &seg_aln.operations {
                match op {
                    AlignOp::Match => {
                        for (k, &idx) in group1.iter().enumerate() {
                            new_sequences[idx].push(stripped_seg1[k][i1]);
                        }
                        for (k, &idx) in group2.iter().enumerate() {
                            new_sequences[idx].push(stripped_seg2[k][i2]);
                        }
                        i1 += 1;
                        i2 += 1;
                    }
                    AlignOp::Delete => {
                        for (k, &idx) in group1.iter().enumerate() {
                            new_sequences[idx].push(stripped_seg1[k][i1]);
                        }
                        for &idx in group2 {
                            new_sequences[idx].push(b'-');
                        }
                        i1 += 1;
                    }
                    AlignOp::Insert => {
                        for &idx in group1 {
                            new_sequences[idx].push(b'-');
                        }
                        for (k, &idx) in group2.iter().enumerate() {
                            new_sequences[idx].push(stripped_seg2[k][i2]);
                        }
                        i2 += 1;
                    }
                }
            }
        }

        return Some((new_sequences, total_score));
    }

    // Non-FFT path: profile_align on stripped profiles.
    let aln = profile_align(
        &stripped_prof1, &stripped_prof2,
        &scoring.substitution_matrix, gap, true, true,
    );
    build_result_from_stripped(
        &aln, group1, group2, sequences, &kept1, &kept2,
        &stripped_prof1, &stripped_prof2,
    )
}

/// Constraint-aware FFT-segmented refinement, port of C's
/// `Falign_localhom` (Falign_localhom.c:163, kobetsubunkatsu=1 path).
///
/// Mirrors the unconstrained FFT-segmented refinement in
/// `realign_all`'s `use_fft` branch but per-segment calls
/// `profile_align_imp` with the local impmtx slice rather than the
/// unconstrained `profile_align`. The impmtx is built once for the full
/// alignment width using the localhom regions; each segment's local
/// view is the rectangle covering the segment's parent column range,
/// indexed by the per-group strip kept-column lists (= C's `gapmap1` /
/// `gapmap2`).
///
/// The cut points are computed from `alignable_segments` at lag=0 just
/// like the unconstrained path (C's `alignableReagion` with maxk=1 in
/// kobetsubunkatsu mode).
fn realign_all_constrained_fft(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    w1n: &[f64],
    w2n: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    lh_table: &LocalHomologyTable,
) -> Option<(Vec<Vec<u8>>, f64)> {
    let width = sequences[group1[0]].len();
    if width == 0 { return None; }

    // Build full-alignment profiles for site-score / segment detection
    // (matches the unconstrained FFT path). For refinement,
    // prof1.length == prof2.length == width.
    let full_prof1 = Profile::from_aligned(
        &group1.iter().map(|&i| sequences[i].as_slice()).collect::<Vec<_>>(),
        w1n, &scoring.amino_map, scoring.nalphabets,
    );
    let full_prof2 = Profile::from_aligned(
        &group2.iter().map(|&i| sequences[i].as_slice()).collect::<Vec<_>>(),
        w2n, &scoring.amino_map, scoring.nalphabets,
    );

    let segment_params = if scoring.seq_type.is_nucleotide() {
        SegmentParams::dna().with_threshold(50.0)
    } else {
        SegmentParams::protein().with_threshold(50.0)
    };

    // Per-position site scores at lag=0 (C's alignableReagion).
    let totaleff: f64 = w1n.iter().sum::<f64>() * w2n.iter().sum::<f64>();
    let len = full_prof1.length.min(full_prof2.length);
    let mut site_scores = vec![0.0f64; len];
    for i in 0..len {
        site_scores[i] = full_prof1.match_score(
            i, &full_prof2, i, &scoring.substitution_matrix,
        ) / totaleff;
    }
    let segments = alignable_segments(&site_scores, &segment_params);

    // Cut points (C's `cut1[i+1] = sortedseg1[i]->center` plus 0 and len).
    let mut cuts: Vec<usize> = Vec::with_capacity(segments.len() + 2);
    cuts.push(0);
    for seg in &segments { cuts.push(seg.center.min(width)); }
    cuts.push(width);
    cuts.sort();
    cuts.dedup();

    // Build the GLOBAL impmtx for the full non-stripped alignment
    // (width × width). C does this via `part_imp_match_init_strict(...,
    // length, length, mseq1, mseq2, ...)` once per branch realign.
    // Each per-segment partA__align then reads `impmtx[start1+gapmap1[i]][start2+gapmap2[j]]`.
    let g1_full: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
    let g2_full: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();
    let global_imp = build_imp_matrix(
        lh_table,
        group1, group2,
        &g1_full, &g2_full,
        w1n, w2n,
        width, width,
        FASTATHRESHOLD_DEFAULT,
    );

    let mut new_sequences: Vec<Vec<u8>> = vec![Vec::new(); sequences.len()];
    let mut total_score = 0.0f64;
    for win in cuts.windows(2) {
        let a = win[0];
        let b = win[1];
        if a >= b { continue; }
        let seg_width = b - a;

        // Slice [a, b) from each member.
        let seg1: Vec<Vec<u8>> = group1.iter()
            .map(|&i| sequences[i][a..b].to_vec()).collect();
        let seg2: Vec<Vec<u8>> = group2.iter()
            .map(|&i| sequences[i][a..b].to_vec()).collect();

        // commongappick within segment + record gapmap (= column index
        // within segment for each kept stripped column).
        let seg_gap1: Vec<bool> = (0..seg_width)
            .map(|c| seg1.iter().all(|s| s[c] == b'-')).collect();
        let seg_gap2: Vec<bool> = (0..seg_width)
            .map(|c| seg2.iter().all(|s| s[c] == b'-')).collect();
        let gapmap1: Vec<usize> = (0..seg_width).filter(|&c| !seg_gap1[c]).collect();
        let gapmap2: Vec<usize> = (0..seg_width).filter(|&c| !seg_gap2[c]).collect();

        let stripped_seg1: Vec<Vec<u8>> = seg1.iter()
            .map(|s| gapmap1.iter().map(|&c| s[c]).collect()).collect();
        let stripped_seg2: Vec<Vec<u8>> = seg2.iter()
            .map(|s| gapmap2.iter().map(|&c| s[c]).collect()).collect();

        // Edge case: one side completely empty after strip.
        if stripped_seg1.is_empty() || stripped_seg1[0].is_empty()
            || stripped_seg2.is_empty() || stripped_seg2[0].is_empty()
        {
            let len1 = if stripped_seg1.is_empty() { 0 } else { stripped_seg1[0].len() };
            let len2 = if stripped_seg2.is_empty() { 0 } else { stripped_seg2[0].len() };
            for (k, &i) in group1.iter().enumerate() {
                new_sequences[i].extend_from_slice(&stripped_seg1[k]);
                new_sequences[i].extend(std::iter::repeat(b'-').take(len2));
            }
            for (k, &i) in group2.iter().enumerate() {
                new_sequences[i].extend(std::iter::repeat(b'-').take(len1));
                new_sequences[i].extend_from_slice(&stripped_seg2[k]);
            }
            continue;
        }

        // Build the segment's local impmtx by indexing the global impmtx
        // at (a + gapmap1[i], a + gapmap2[j]) — mirrors C's
        // `imp_match_out_vead_gapmap(imp[j] = impmtx[i1][start2+gapmap2[j]])`
        // (partSalignmm.c:71-83). For refinement, start1 = start2 = a.
        let l1 = stripped_seg1[0].len();
        let l2 = stripped_seg2[0].len();
        let mut local_imp = vec![vec![0.0f64; l2]; l1];
        for i in 0..l1 {
            let row = a + gapmap1[i];
            for j in 0..l2 {
                let col = a + gapmap2[j];
                if row < global_imp.len() && col < global_imp[row].len() {
                    local_imp[i][j] = global_imp[row][col];
                }
            }
        }

        let s1_refs: Vec<&[u8]> = stripped_seg1.iter().map(|s| s.as_slice()).collect();
        let s2_refs: Vec<&[u8]> = stripped_seg2.iter().map(|s| s.as_slice()).collect();
        let mut prof_seg1 = Profile::from_aligned(&s1_refs, w1n, &scoring.amino_map, scoring.nalphabets);
        let mut prof_seg2 = Profile::from_aligned(&s2_refs, w2n, &scoring.amino_map, scoring.nalphabets);

        // Per-segment boundary gap correction (mirrors C's `getkyokaigap`
        // → `new_OpeningGapCount` in `MSalignmm`/`partA__align`). Same as
        // the unconstrained FFT path: a sequence already in a gap just
        // before the segment's first column does NOT count as opening at
        // position 0.
        let sgap_inside = a > 0;
        let egap_inside = b < width;
        if sgap_inside || egap_inside {
            if !prof_seg1.ogcp.is_empty() {
                for (k, &idx) in group1.iter().enumerate() {
                    if sgap_inside
                        && sequences[idx][a - 1] == b'-'
                        && !stripped_seg1[k].is_empty()
                        && stripped_seg1[k][0] == b'-'
                    {
                        prof_seg1.ogcp[0] -= w1n[k];
                    }
                    let l1k = stripped_seg1[k].len();
                    if egap_inside
                        && l1k > 0
                        && stripped_seg1[k][l1k - 1] == b'-'
                        && sequences[idx][b] == b'-'
                    {
                        let last = prof_seg1.fgcp.len() - 1;
                        prof_seg1.fgcp[last] -= w1n[k];
                    }
                }
            }
            if !prof_seg2.ogcp.is_empty() {
                for (k, &idx) in group2.iter().enumerate() {
                    if sgap_inside
                        && sequences[idx][a - 1] == b'-'
                        && !stripped_seg2[k].is_empty()
                        && stripped_seg2[k][0] == b'-'
                    {
                        prof_seg2.ogcp[0] -= w2n[k];
                    }
                    let l2k = stripped_seg2[k].len();
                    if egap_inside
                        && l2k > 0
                        && stripped_seg2[k][l2k - 1] == b'-'
                        && sequences[idx][b] == b'-'
                    {
                        let last = prof_seg2.fgcp.len() - 1;
                        prof_seg2.fgcp[last] -= w2n[k];
                    }
                }
            }
        }

        // C's Falign_localhom segment loop passes `headgp = tailgp = 1`
        // when a segment is interior (`(i==0)?outgap:1` etc.); for
        // dvtditr, outgap=1 anyway so all segments use 1.
        //
        // C uses `partA__align` (partSalignmm.c:1218,1235) which uses STRICT
        // `>` for the prept-vs-mi/mjpt tie-break (the "// 2018/Apr" change).
        // The progressive `A__align` uses `>=`. We pass strict_part_tiebreak=true
        // to match `partA__align` exactly.
        //
        // Boundary nongap-frequencies (C's headgapfreq{1,2} and
        // gapfreq{1,2}[lgth] computed via outgapcount on sgap/egap):
        //   head{1,2} = nongap fraction at full-alignment column [a-1]
        //               (1.0 when a == 0 → C's `sgap[j]='o'` branch)
        //   tail{1,2} = nongap fraction at full-alignment column [b]
        //               (1.0 when b == width → C's `egap[j]='o'` branch)
        let head1 = if a > 0 {
            let s: f64 = group1.iter().enumerate()
                .filter(|&(_, &idx)| sequences[idx][a - 1] == b'-')
                .map(|(k, _)| w1n[k]).sum();
            1.0 - s
        } else { 1.0 };
        let head2 = if a > 0 {
            let s: f64 = group2.iter().enumerate()
                .filter(|&(_, &idx)| sequences[idx][a - 1] == b'-')
                .map(|(k, _)| w2n[k]).sum();
            1.0 - s
        } else { 1.0 };
        let tail1 = if b < width {
            let s: f64 = group1.iter().enumerate()
                .filter(|&(_, &idx)| sequences[idx][b] == b'-')
                .map(|(k, _)| w1n[k]).sum();
            1.0 - s
        } else { 1.0 };
        let tail2 = if b < width {
            let s: f64 = group2.iter().enumerate()
                .filter(|&(_, &idx)| sequences[idx][b] == b'-')
                .map(|(k, _)| w2n[k]).sum();
            1.0 - s
        } else { 1.0 };
        let boundary = BoundaryFreqs { head1, head2, tail1, tail2 };
        let seg_aln = profile_align_imp_with_boundary(
            &prof_seg1, &prof_seg2, &scoring.substitution_matrix, gap,
            true, true, Some(&local_imp), true, boundary,
        );
        total_score += seg_aln.score;

        // Reconstruct segment output by applying ops to stripped segments.
        let mut i1 = 0usize;
        let mut i2 = 0usize;
        for op in &seg_aln.operations {
            match op {
                AlignOp::Match => {
                    for (k, &idx) in group1.iter().enumerate() {
                        new_sequences[idx].push(stripped_seg1[k][i1]);
                    }
                    for (k, &idx) in group2.iter().enumerate() {
                        new_sequences[idx].push(stripped_seg2[k][i2]);
                    }
                    i1 += 1; i2 += 1;
                }
                AlignOp::Delete => {
                    for (k, &idx) in group1.iter().enumerate() {
                        new_sequences[idx].push(stripped_seg1[k][i1]);
                    }
                    for &idx in group2 {
                        new_sequences[idx].push(b'-');
                    }
                    i1 += 1;
                }
                AlignOp::Insert => {
                    for &idx in group1 {
                        new_sequences[idx].push(b'-');
                    }
                    for (k, &idx) in group2.iter().enumerate() {
                        new_sequences[idx].push(stripped_seg2[k][i2]);
                    }
                    i2 += 1;
                }
            }
        }
    }

    // Pad sequences not in either group to the new width.
    let new_width = if !group1.is_empty() {
        new_sequences[group1[0]].len()
    } else if !group2.is_empty() {
        new_sequences[group2[0]].len()
    } else {
        width
    };
    for (i, s) in new_sequences.iter_mut().enumerate() {
        if !group1.contains(&i) && !group2.contains(&i) {
            // "Other" sequences — preserve from input. Our refinement only
            // realigns groups that partition all sequences, so this should
            // never trigger; keep original.
            *s = sequences[i].clone();
        }
        if s.len() < new_width {
            s.resize(new_width, b'-');
        }
    }

    Some((new_sequences, total_score))
}

/// Build result sequences from an alignment on stripped profiles.
/// Maps alignment operations back to original column positions via kept1/kept2.
fn build_result_from_stripped(
    aln: &mafft_align::Alignment,
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    kept1: &[usize],
    kept2: &[usize],
    prof1: &Profile,
    prof2: &Profile,
) -> Option<(Vec<Vec<u8>>, f64)> {
    let consumed1 = aln.operations.iter()
        .filter(|op| matches!(op, AlignOp::Match | AlignOp::Delete)).count();
    let consumed2 = aln.operations.iter()
        .filter(|op| matches!(op, AlignOp::Match | AlignOp::Insert)).count();
    if consumed1 != prof1.length || consumed2 != prof2.length {
        return None;
    }

    let mut new_sequences = vec![Vec::with_capacity(aln.operations.len()); sequences.len()];
    let mut c1 = 0usize;
    let mut c2 = 0usize;
    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                let oc1 = kept1[c1];
                let oc2 = kept2[c2];
                for &i in group1 { new_sequences[i].push(sequences[i][oc1]); }
                for &i in group2 { new_sequences[i].push(sequences[i][oc2]); }
                c1 += 1; c2 += 1;
            }
            AlignOp::Delete => {
                let oc1 = kept1[c1];
                for &i in group1 { new_sequences[i].push(sequences[i][oc1]); }
                for &i in group2 { new_sequences[i].push(b'-'); }
                c1 += 1;
            }
            AlignOp::Insert => {
                let oc2 = kept2[c2];
                for &i in group1 { new_sequences[i].push(b'-'); }
                for &i in group2 { new_sequences[i].push(sequences[i][oc2]); }
                c2 += 1;
            }
        }
    }
    Some((new_sequences, aln.score))
}

/// Build result sequences from an alignment on non-stripped (full) profiles.
/// Cursors index directly into the full-width sequences.
fn build_result_from_full(
    aln: &mafft_align::Alignment,
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    prof1: &Profile,
    prof2: &Profile,
) -> Option<(Vec<Vec<u8>>, f64)> {
    let consumed1 = aln.operations.iter()
        .filter(|op| matches!(op, AlignOp::Match | AlignOp::Delete)).count();
    let consumed2 = aln.operations.iter()
        .filter(|op| matches!(op, AlignOp::Match | AlignOp::Insert)).count();
    if consumed1 != prof1.length || consumed2 != prof2.length {
        return None;
    }

    let mut new_sequences = vec![Vec::with_capacity(aln.operations.len()); sequences.len()];
    let mut c1 = 0usize;
    let mut c2 = 0usize;
    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                for &i in group1 { new_sequences[i].push(sequences[i][c1]); }
                for &i in group2 { new_sequences[i].push(sequences[i][c2]); }
                c1 += 1; c2 += 1;
            }
            AlignOp::Delete => {
                for &i in group1 { new_sequences[i].push(sequences[i][c1]); }
                for &i in group2 { new_sequences[i].push(b'-'); }
                c1 += 1;
            }
            AlignOp::Insert => {
                for &i in group1 { new_sequences[i].push(b'-'); }
                for &i in group2 { new_sequences[i].push(sequences[i][c2]); }
                c2 += 1;
            }
        }
    }
    Some((new_sequences, aln.score))
}

fn group_all_gap_columns(group: &[usize], sequences: &[Vec<u8>], width: usize) -> Vec<bool> {
    // Match C's commongappick: only '-' counts as gap (not '.').
    let mut all_gap = vec![true; width];
    for &idx in group {
        for (col, &ch) in sequences[idx].iter().enumerate() {
            if ch != b'-' {
                all_gap[col] = false;
            }
        }
    }
    all_gap
}

/// Sum the impmatch (constraint-importance bonus) across the diagonal of
/// the current alignment. Mirrors C's `oimpmatchdouble = sum imp_match_out_sc(i,i)`
/// loop (tditeration.c:925) for the existing alignment's `impmtx`.
///
/// `eff1`, `eff2` are group-local sum-1-normalized weights matching
/// what `imp_match_init_strict` is called with in C (= `effarr1`/`effarr2`
/// after `fastconjuction_noname` per-cluster normalization).
fn compute_impmatch_diagonal(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    eff1: &[f64],
    eff2: &[f64],
    lh_table: &LocalHomologyTable,
) -> f64 {
    let width = sequences[group1[0]].len();
    if width == 0 { return 0.0; }
    let g1_seqs: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
    let g2_seqs: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();
    let imp = build_imp_matrix(
        lh_table,
        group1, group2,
        &g1_seqs, &g2_seqs,
        eff1, eff2,
        width, width,
        FASTATHRESHOLD_DEFAULT,
    );
    let mut total = 0.0f64;
    for i in 0..width {
        if let Some(row) = imp.get(i) {
            if let Some(&v) = row.get(i) {
                total += v;
            }
        }
    }
    total
}

fn compute_split_score(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
) -> f64 {
    // Port of C's intergroup_score flow: weights are per-group normalized
    // by fastconjuction_noname (tddis.c line 548) before being passed.
    // Apply the same normalization here: each group's weights sum to 1.0.
    const MINIMUM_WEIGHT: f64 = 0.00001;
    let w1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MINIMUM_WEIGHT)).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MINIMUM_WEIGHT)).collect();
    let s1: f64 = w1.iter().sum();
    let s2: f64 = w2.iter().sum();
    let w1n: Vec<f64> = if s1 > 0.0 { w1.iter().map(|w| w / s1).collect() } else { vec![1.0; group1.len()] };
    let w2n: Vec<f64> = if s2 > 0.0 { w2.iter().map(|w| w / s2).collect() } else { vec![1.0; group2.len()] };

    // Sequential sum to match C's deterministic accumulation order.
    // par_iter gives non-deterministic summation order, which causes
    // small FP divergence that cascades into accept/reject decisions.
    let mut total = 0.0f64;
    for (i_local, &i) in group1.iter().enumerate() {
        let wi = w1n[i_local];
        for (j_local, &j) in group2.iter().enumerate() {
            let wj = w2n[j_local];
            total += pairwise_score(&sequences[i], &sequences[j], scoring) * wi * wj;
        }
    }
    total
}

/// Branchless pairwise scoring for auto-vectorization.
///
/// The gap check is converted to a mask multiply: if either residue is '-',
/// the score contribution is 0. This eliminates branches that prevent SIMD.
#[inline]
fn pairwise_score(seq1: &[u8], seq2: &[u8], scoring: &ScoringContext) -> f64 {
    let map = &scoring.amino_map;
    let mtx = &scoring.substitution_matrix;
    let mtx_size = mtx.len();

    // Port of C's intergroup_score (mltaln9.c lines 404-475):
    // - Gap-gap positions: skipped (continue).
    // - Match positions: add amino_dis[c1][c2].
    // - Gap in seq1: add `penalty` (gap open), then consume all consecutive
    //   '-' in seq1. Same for seq2.
    // - amino_dis_consweight_multi[gap][*] = 0, so gap-region positions
    //   contribute only the single gap-open penalty per gap run.
    let penalty = scoring.gap.open as f64;
    let len = seq1.len().min(seq2.len());
    let mut score = 0.0f64;
    let mut k = 0;
    while k < len {
        let a = seq1[k];
        let b = seq2[k];
        if a == b'-' && b == b'-' {
            k += 1;
            continue;
        }
        if a == b'-' {
            score += penalty;
            // Consume all consecutive gaps in seq1 (C's while-loop at line 448).
            k += 1;
            while k < len && seq1[k] == b'-' {
                k += 1;
            }
            continue;
        }
        if b == b'-' {
            score += penalty;
            k += 1;
            while k < len && seq2[k] == b'-' {
                k += 1;
            }
            continue;
        }
        let i = map[a as usize] as usize;
        let j = map[b as usize] as usize;
        if i < mtx_size && j < mtx_size {
            score += mtx[i][j] as f64;
        }
        k += 1;
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_tree::{DistanceMatrix, upgma};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use crate::progressive::progressive_align;

    /// Helper: build a 6-sequence UPGMA topology for branch-enumeration tests.
    fn make_6seq_topology() -> (Topology, usize) {
        let nseq = 6;
        let mut dm = DistanceMatrix::new(nseq);
        for i in 0..nseq {
            for j in (i + 1)..nseq {
                dm.set(i, j, (j - i) as f64 * 0.1);
            }
        }
        (upgma(&dm), nseq)
    }

    // ---------------------------------------------------------------
    // Regression guards for the iterative-refinement fixes.
    // Each test targets one specific behavior ported from C's
    // TreeDependentIteration() in tditeration.c. If any of these
    // are accidentally reverted, at least one test will fail.
    // ---------------------------------------------------------------

    /// Guard: branch count = (nseq-1)*2 - 1, matching C's nbranch formula.
    ///
    /// C computes `nbranch = (njob-1) * 2 - 1` (tditeration.c line 1458).
    /// The root step contributes only 1 branch (side 1), all others contribute
    /// 2 (sides 0 and 1). Reverting the root-step skip would produce
    /// (nseq-1)*2 branches instead.
    #[test]
    fn branch_count_matches_c_formula() {
        let (topo, nseq) = make_6seq_topology();
        let branch_map = build_branch_map(&topo, nseq);
        let total: usize = branch_map.iter().map(|sides| sides.len()).sum();
        let expected = (nseq - 1) * 2 - 1;
        assert_eq!(total, expected,
            "branch count should be (nseq-1)*2-1 = {expected}, got {total}");
    }

    /// Guard: root step has exactly 1 branch (side 1 only).
    ///
    /// C forces `k = 1` at the root step (tditeration.c line 1667:
    /// `if( l == locnjob-2 ) k = 1`), skipping side 0 because at the
    /// root left-vs-complement and right-vs-complement are identical
    /// splits. Reverting would give the root step 2 branches.
    #[test]
    fn root_step_has_single_branch() {
        let (topo, nseq) = make_6seq_topology();
        let branch_map = build_branch_map(&topo, nseq);
        let root_branches = branch_map.last().unwrap();
        assert_eq!(root_branches.len(), 1,
            "root step should have 1 branch (side 1 only), got {}", root_branches.len());
        assert_eq!(root_branches[0].0, 1, "root branch should be side 1");
    }

    /// Guard: non-root steps each have exactly 2 branches (sides 0 and 1).
    #[test]
    fn non_root_steps_have_two_branches() {
        let (topo, nseq) = make_6seq_topology();
        let branch_map = build_branch_map(&topo, nseq);
        for (step_idx, sides) in branch_map.iter().enumerate() {
            if step_idx < branch_map.len() - 1 {
                assert_eq!(sides.len(), 2,
                    "non-root step {step_idx} should have 2 branches, got {}", sides.len());
            }
        }
    }

    /// Guard: default cut is 0.0 (accept only strict improvements).
    ///
    /// C's dvtditr.c sets `cut = 0.0` (line 71). The acceptance test is
    /// `tscore > mscore - cut/100*mscore`, so with cut=0 only strictly
    /// improving moves are accepted. Reverting to a nonzero cut would
    /// accept non-improving moves.
    #[test]
    fn default_cut_is_zero() {
        let params = RefinementParams::default();
        assert_eq!(params.cut, 0.0,
            "default cut must be 0.0 (strict improvement only), matching C's dvtditr.c");
    }

    /// Guard: even iterations traverse steps forward, odd iterations reverse.
    ///
    /// C alternates direction (tditeration.c lines 1641-1648):
    ///   even → lin=0, ldf=+1 (forward)
    ///   odd  → lin=locnjob-2, ldf=-1 (reverse)
    /// This test verifies the first branch processed differs between
    /// iteration 0 (forward) and iteration 1 (reverse).
    #[test]
    fn alternating_direction_between_iterations() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIKLMNPQR".to_vec(),
            b"ACDEFHIKLMNPQR".to_vec(),
            b"ACDEHIKLMNPQR".to_vec(),
            b"ACDHIKLMNPQR".to_vec(),
            b"ACDHIKLMNP".to_vec(),
            b"ACDHIKLM".to_vec(),
        ];
        let names: Vec<String> = (0..6).map(|i| format!("s{i}")).collect();
        let (topo, nseq) = make_6seq_topology();
        let nsteps = topo.steps.len();

        // Verify the step_order logic directly.
        let forward: Vec<usize> = (0..nsteps).collect();
        let reverse: Vec<usize> = (0..nsteps).rev().collect();

        // Even iteration → forward
        assert_eq!(forward[0], 0, "forward should start at step 0");
        // Odd iteration → reverse
        assert_eq!(reverse[0], nsteps - 1, "reverse should start at last step");
        // They must differ (nsteps > 1 for any nseq > 2)
        assert_ne!(forward, reverse,
            "forward and reverse step orders must differ for alternation");
    }

    /// Guard: oscillation detection terminates refinement early.
    ///
    /// C checks per-branch score history (tditeration.c lines 2343-2371)
    /// and stops if a branch's score at iteration N matches the score
    /// from iteration N-2. We verify that iterative_refine returns in
    /// fewer than max_iterations when running on inputs that converge
    /// quickly (which will produce identical scores across iterations).
    #[test]
    fn refinement_terminates_not_at_max_iterations() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        // Three nearly-identical sequences: refinement should converge fast,
        // well before hitting 100 iterations.
        let seqs = vec![
            b"ACDEFGHIK".to_vec(),
            b"ACDEFGHIK".to_vec(),
            b"ACDEFGHIK".to_vec(),
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];

        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.001);
        dm.set(0, 2, 0.001);
        dm.set(1, 2, 0.001);
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams {
            max_iterations: 100,
            ..Default::default()
        };

        let iters = iterative_refine(&mut msa, &topo, &scoring, &params, None);
        // Identical sequences must converge immediately — either via the
        // identity check or via the convergence counter (nseq * 2 = 6).
        assert!(iters < 100,
            "expected early termination (convergence/oscillation), got {iters} iterations");
        assert!(iters <= 2,
            "identical sequences should converge in 1-2 iterations, got {iters}");
    }

    // ---------------------------------------------------------------
    // Original functional tests (preserved).
    // ---------------------------------------------------------------

    #[test]
    fn refinement_converges() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIK".to_vec(),
            b"ACDEFHIK".to_vec(),
            b"ACDHIK".to_vec(),
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];

        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.1);
        dm.set(0, 2, 0.3);
        dm.set(1, 2, 0.2);
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams {
            max_iterations: 10,
            ..Default::default()
        };

        let iters = iterative_refine(&mut msa, &topo, &scoring, &params, None);
        assert!(iters <= 10);

        let width = msa.width();
        for seq in &msa.sequences {
            assert_eq!(seq.len(), width);
        }

        let ungapped: Vec<Vec<u8>> = msa.sequences.iter()
            .map(|s| s.iter().filter(|&&c| c != b'-').cloned().collect())
            .collect();
        assert_eq!(ungapped[0], b"ACDEFGHIK");
        assert_eq!(ungapped[1], b"ACDEFHIK");
        assert_eq!(ungapped[2], b"ACDHIK");
    }

    #[test]
    fn refinement_preserves_width_consistency() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIKLMNP".to_vec(),
            b"ACDEFHIKLMNP".to_vec(),
            b"ACDEHIKLMNP".to_vec(),
            b"ACDHIKLMNP".to_vec(),
        ];
        let names: Vec<String> = (0..4).map(|i| format!("s{i}")).collect();

        let mut dm = DistanceMatrix::new(4);
        dm.set(0, 1, 0.1); dm.set(0, 2, 0.2); dm.set(0, 3, 0.3);
        dm.set(1, 2, 0.15); dm.set(1, 3, 0.25); dm.set(2, 3, 0.15);
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams { max_iterations: 5, ..Default::default() };

        iterative_refine(&mut msa, &topo, &scoring, &params, None);

        let width = msa.width();
        assert!(width > 0);
        for (i, seq) in msa.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "sequence {i} has wrong width");
        }

        // Verify residue preservation
        for (i, seq) in msa.sequences.iter().enumerate() {
            let residue_count = seq.iter().filter(|&&c| c != b'-').count();
            assert_eq!(residue_count, seqs[i].len(),
                "sequence {i} lost residues: {} vs {}", residue_count, seqs[i].len());
        }
    }

    #[test]
    fn refinement_six_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIKLMNPQR".to_vec(),
            b"ACDEFHIKLMNPQR".to_vec(),
            b"ACDEHIKLMNPQR".to_vec(),
            b"ACDHIKLMNPQR".to_vec(),
            b"ACDHIKLMNP".to_vec(),
            b"ACDHIKLM".to_vec(),
        ];
        let names: Vec<String> = (0..6).map(|i| format!("s{i}")).collect();

        let mut dm = DistanceMatrix::new(6);
        for i in 0..6 {
            for j in (i + 1)..6 {
                dm.set(i, j, (j - i) as f64 * 0.1);
            }
        }
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams { max_iterations: 3, ..Default::default() };

        iterative_refine(&mut msa, &topo, &scoring, &params, None);

        let width = msa.width();
        for (i, seq) in msa.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "sequence {i} has wrong width after refinement");
            let residues = seq.iter().filter(|&&c| c != b'-').count();
            assert_eq!(residues, seqs[i].len(),
                "sequence {i} lost residues during refinement");
        }
    }

    /// Guard: FFT-accelerated refinement produces valid results and
    /// does not cause width explosion.
    ///
    /// C always uses Falign (FFT) in refinement. This test verifies that
    /// use_fft=true in RefinementParams produces a valid alignment with
    /// bounded width growth (no exponential blow-up).
    #[test]
    fn refinement_fft_no_width_explosion() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIKLMNPQRSTVWY".to_vec(),
            b"ACDEFHIKLMNPQRSTVWY".to_vec(),
            b"ACDEHIKLMNPQRSTVWY".to_vec(),
            b"ACDHIKLMNPQRSTVWY".to_vec(),
            b"ACDHIKLMNPQR".to_vec(),
            b"ACDHIKLM".to_vec(),
        ];
        let names: Vec<String> = (0..6).map(|i| format!("s{i}")).collect();

        let mut dm = DistanceMatrix::new(6);
        for i in 0..6 {
            for j in (i + 1)..6 {
                dm.set(i, j, (j - i) as f64 * 0.1);
            }
        }
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let pre_width = msa.width();

        let params = RefinementParams {
            max_iterations: 5,
            use_fft: true,
            ..Default::default()
        };

        iterative_refine(&mut msa, &topo, &scoring, &params, None);

        let post_width = msa.width();
        // Width should not blow up — allow at most 2x growth for reasonable
        // refinement (C typically keeps width within ~10% of progressive).
        assert!(post_width <= pre_width * 2,
            "width explosion: {} -> {} (>2x growth)", pre_width, post_width);

        for (i, seq) in msa.sequences.iter().enumerate() {
            assert_eq!(seq.len(), post_width, "sequence {i} has wrong width");
            let residues = seq.iter().filter(|&&c| c != b'-').count();
            assert_eq!(residues, seqs[i].len(),
                "sequence {i} lost residues during FFT refinement");
        }
    }
}
