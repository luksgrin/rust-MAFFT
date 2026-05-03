/// Progressive alignment following a guide tree.
///
/// Matches C's `treebase()` from `disttbfast.c`: at each merge step,
/// only group1 and group2 sequences are modified. "Other" sequences
/// are left untouched. Profiles are cached after each merge step
/// (`cpmxhist`) to match C's exact float accumulation order.

use std::collections::HashMap;
use mafft_align::{profile_align, pairwise_align11, fft_profile_align, Profile, GapModel, Alignment, AlignOp, FftAlignParams};
use mafft_tree::{Topology, sequence_weights};
use mafft_types::ScoringContext;

#[derive(Debug, Clone)]
pub struct MultipleAlignment {
    pub sequences: Vec<Vec<u8>>,
    pub names: Vec<String>,
    pub score: f64,
    /// Per-merge-step trace: one entry per progressive merge. Each entry is
    /// `(clus1_size, clus2_size, width_after_merge, score)` matching C's
    /// `RDBG step clus1 clus2 width score` debug line. Used for regression
    /// tests that assert byte-level parity with C on a per-step basis;
    /// harmless to ignore.
    pub step_trace: Vec<StepTrace>,
}

#[derive(Debug, Clone, Copy)]
pub struct StepTrace {
    pub clus1: usize,
    pub clus2: usize,
    pub width: usize,
    pub score: f64,
}

impl MultipleAlignment {
    pub fn width(&self) -> usize {
        self.sequences.first().map_or(0, |s| s.len())
    }
    pub fn nseq(&self) -> usize {
        self.sequences.len()
    }
}

/// Cached profile from a previous merge step.
/// Stores the blended composition probability matrix, gap frequencies,
/// and opening/closing gap counts — matching C's `cpmxhist`.
#[derive(Debug, Clone)]
struct CachedProfile {
    profile: Profile,
    /// Effective weight of this group (orieff), for blending at the next merge.
    eff: f64,
}

pub fn progressive_align(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
) -> MultipleAlignment {
    progressive_align_with_constraints(
        sequences, names, topology, scoring, use_fft, shift_penalty, None,
    )
}

/// Progressive alignment with optional local-homology constraints.
///
/// When `constraints` is `Some`, every merge calls `profile_align_imp`
/// with a per-merge impmtx built via `mafft-align::build_imp_matrix` —
/// mirroring what `tbfast` does in C's L-INS-i pipeline (`Falign_localhom`,
/// or the non-FFT `partA__align` per segment). Without this the initial
/// progressive alignment is FFT-NS-i-like and only refinement sees the
/// constraints, leaving L-INS-i/E-INS-i shapes systematically off vs C.
pub fn progressive_align_with_constraints(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
) -> MultipleAlignment {
    let nseq = sequences.len();
    if nseq == 0 {
        return MultipleAlignment {
            sequences: Vec::new(), names: Vec::new(), score: 0.0, step_trace: Vec::new(),
        };
    }
    if nseq == 1 {
        return MultipleAlignment {
            sequences: sequences.to_vec(), names: names.to_vec(), score: 0.0, step_trace: Vec::new(),
        };
    }

    let weights = sequence_weights(topology);
    let mut aligned: Vec<Vec<u8>> = sequences.to_vec();

    let mut last_score = 0.0;
    let mut gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    if let Some(shift) = shift_penalty {
        gap = gap.with_shift(shift);
    }

    // Profile cache: maps a set of sequence indices (sorted) to its cached profile.
    // After each merge, the merged profile is stored so the next merge can reuse it.
    let mut profile_cache: HashMap<Vec<usize>, CachedProfile> = HashMap::new();

    let mut step_trace: Vec<StepTrace> = Vec::with_capacity(topology.steps.len());
    for (step_idx, step) in topology.steps.iter().enumerate() {
        last_score = merge_step_cached(
            &step.left, &step.right, &mut aligned, &weights, scoring, &gap, use_fft,
            &mut profile_cache, constraints,
        );

        let width = aligned[step.left[0]].len().max(aligned[step.right[0]].len());
        step_trace.push(StepTrace {
            clus1: step.left.len(),
            clus2: step.right.len(),
            width,
            score: last_score,
        });
        if std::env::var("MAFFT_DEBUG_STEPS").is_ok() {
            eprintln!("RDBG {} {} {} {} {:.1}",
                step_idx, step.left.len(), step.right.len(), width, last_score);
        }
    }

    let max_width = aligned.iter().map(|s| s.len()).max().unwrap_or(0);
    for seq in &mut aligned {
        seq.resize(max_width, b'-');
    }

    MultipleAlignment {
        sequences: aligned, names: names.to_vec(), score: last_score, step_trace,
    }
}

fn merge_step_cached(
    group1: &[usize],
    group2: &[usize],
    aligned: &mut Vec<Vec<u8>>,
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    use_fft: bool,
    cache: &mut HashMap<Vec<usize>, CachedProfile>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
) -> f64 {
    let width1 = aligned[group1[0]].len();
    let width2 = aligned[group2[0]].len();

    // Look up cached profiles or build from sequences
    let key1 = sorted_key(group1);
    let key2 = sorted_key(group2);

    let (prof1, eff1) = if let Some(cached) = cache.get(&key1) {
        (cached.profile.clone(), cached.eff)
    } else {
        let (prof, eff) = build_profile_from_seqs(group1, aligned, weights, scoring);
        (prof, eff)
    };

    let (prof2, eff2) = if let Some(cached) = cache.get(&key2) {
        (cached.profile.clone(), cached.eff)
    } else {
        let (prof, eff) = build_profile_from_seqs(group2, aligned, weights, scoring);
        (prof, eff)
    };

    // C uses Falign (FFT-accelerated) for ALL steps when ffttry is true
    // (nlen > clus, which is always true). G__align11 is only used when
    // FFT is disabled (use_fft=false) and both groups are single sequences.
    // When alg='A', the non-FFT fallback is A__align (= profile_align).
    //
    // C's disttbfast passes `outgap, outgap` for headgp/tailgp in G__align11
    // and A__align. With the -O flag (always set by mafft script), outgap=0,
    // which means no penalty is applied to terminal gaps (TERMGAPFAC=0).
    let aln = if !use_fft && group1.len() == 1 && group2.len() == 1
        && constraints.is_none()
    {
        // G__align11 path: flat gap penalty, character-level scoring
        // Only used when FFT is disabled and no constraints. With constraints
        // (L-INS-i / E-INS-i / G-INS-i), single-vs-single merges still need
        // the impmtx contribution per cell — fall through to the constrained
        // profile DP below.
        pairwise_align11(
            &aligned[group1[0]], &aligned[group2[0]],
            &scoring.substitution_matrix, &scoring.amino_map,
            scoring.gap.open as f64, false, false,
        )
    } else if use_fft {
        // C uses Falign for ALL steps when use_fft=true (ffttry = nlen > clus,
        // always true). No minimum profile length check.
        let fft_params = FftAlignParams {
            num_candidates: 20,
            segment_params: if scoring.seq_type.is_nucleotide() {
                mafft_fft::SegmentParams::dna()
            } else {
                mafft_fft::SegmentParams::protein()
            },
            gap: gap.clone(),
            head_gap: false,
            tail_gap: false,
            num_channels: scoring.nscoredalphabets,
        };
        fft_profile_align(&prof1, &prof2, &scoring.substitution_matrix, &fft_params)
    } else if let Some(table) = constraints {
        // Constraint-aware progressive merge (L-INS-i / E-INS-i tbfast path).
        // Build per-cell impmtx from the localhom table over the group split,
        // then call the importance-aware DP. Mirrors C's `partA__align` /
        // `Falign_localhom` per-segment DP with `imp_match_out_vead` adding
        // the importance bonus row-by-row.
        let g1_seq_refs: Vec<&[u8]> = group1.iter().map(|&i| aligned[i].as_slice()).collect();
        let g2_seq_refs: Vec<&[u8]> = group2.iter().map(|&i| aligned[i].as_slice()).collect();
        // Group-local sum-1 normalized weights (matches C's
        // `fastconjuction_noname` `peff[m] /= total`, tddis.c:552-556).
        const MINIMUM_WEIGHT: f64 = 0.00001;
        let w1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MINIMUM_WEIGHT)).collect();
        let w2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MINIMUM_WEIGHT)).collect();
        let s1: f64 = w1.iter().sum();
        let s2: f64 = w2.iter().sum();
        let w1n: Vec<f64> = if s1 > 0.0 { w1.iter().map(|w| w / s1).collect() } else { vec![1.0; group1.len()] };
        let w2n: Vec<f64> = if s2 > 0.0 { w2.iter().map(|w| w / s2).collect() } else { vec![1.0; group2.len()] };
        let imp = mafft_align::build_imp_matrix(
            table,
            group1, group2,
            &g1_seq_refs, &g2_seq_refs,
            &w1n, &w2n,
            prof1.length, prof2.length,
            mafft_align::FASTATHRESHOLD_DEFAULT,
        );
        mafft_align::profile_align_imp(
            &prof1, &prof2, &scoring.substitution_matrix, gap,
            false, false, Some(&imp),
        )
    } else {
        profile_align(&prof1, &prof2, &scoring.substitution_matrix, gap, false, false)
    };

    // Build gaptables for profile caching (matching C's gaptable1/gaptable2)
    let new_width = aln.operations.len();
    let mut gaptable1 = Vec::with_capacity(new_width); // 'o' = content, '-' = gap
    let mut gaptable2 = Vec::with_capacity(new_width);
    for op in &aln.operations {
        match op {
            AlignOp::Match => { gaptable1.push(b'o'); gaptable2.push(b'o'); }
            AlignOp::Delete => { gaptable1.push(b'o'); gaptable2.push(b'-'); }
            AlignOp::Insert => { gaptable1.push(b'-'); gaptable2.push(b'o'); }
        }
    }

    // Cache the merged profile (C's createcpmxresult + creategapfreqresult +
    // createogresult + createfgresult). Only cache for groups > 20 sequences
    // (matching C's condition at MSalignmm.c line 2431).
    let total_eff = eff1 + eff2;
    let combined_seqs = group1.len() + group2.len();
    if total_eff > 0.0 && combined_seqs > 20 {
        let norm_eff1 = eff1 / total_eff;
        let norm_eff2 = eff2 / total_eff;
        let merged_prof = blend_profiles_exact(
            &prof1, &prof2,
            norm_eff1, norm_eff2,
            &gaptable1, &gaptable2,
            scoring.nalphabets,
        );
        let mut merged_key = group1.to_vec();
        merged_key.extend_from_slice(group2);
        merged_key.sort();
        cache.insert(merged_key, CachedProfile {
            profile: merged_prof,
            eff: total_eff,
        });
    }

    // Remove child caches (they won't be needed again)
    cache.remove(&key1);
    cache.remove(&key2);

    // Build new sequences for group1 and group2 ONLY
    let mut cursor1 = 0usize;
    let mut cursor2 = 0usize;
    let mut new_seqs_g1: Vec<Vec<u8>> = vec![Vec::with_capacity(new_width); group1.len()];
    let mut new_seqs_g2: Vec<Vec<u8>> = vec![Vec::with_capacity(new_width); group2.len()];

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                for (gi, &idx) in group1.iter().enumerate() {
                    new_seqs_g1[gi].push(if cursor1 < width1 { aligned[idx][cursor1] } else { b'-' });
                }
                for (gi, &idx) in group2.iter().enumerate() {
                    new_seqs_g2[gi].push(if cursor2 < width2 { aligned[idx][cursor2] } else { b'-' });
                }
                cursor1 += 1;
                cursor2 += 1;
            }
            AlignOp::Delete => {
                for (gi, &idx) in group1.iter().enumerate() {
                    new_seqs_g1[gi].push(if cursor1 < width1 { aligned[idx][cursor1] } else { b'-' });
                }
                for gi in 0..group2.len() { new_seqs_g2[gi].push(b'-'); }
                cursor1 += 1;
            }
            AlignOp::Insert => {
                for gi in 0..group1.len() { new_seqs_g1[gi].push(b'-'); }
                for (gi, &idx) in group2.iter().enumerate() {
                    new_seqs_g2[gi].push(if cursor2 < width2 { aligned[idx][cursor2] } else { b'-' });
                }
                cursor2 += 1;
            }
        }
    }

    for (gi, &idx) in group1.iter().enumerate() { aligned[idx] = new_seqs_g1[gi].clone(); }
    for (gi, &idx) in group2.iter().enumerate() { aligned[idx] = new_seqs_g2[gi].clone(); }

    aln.score
}

fn sorted_key(group: &[usize]) -> Vec<usize> {
    let mut k = group.to_vec();
    k.sort();
    k
}

fn build_profile_from_seqs(
    group: &[usize],
    aligned: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
) -> (Profile, f64) {
    let seqs: Vec<&[u8]> = group.iter().map(|&i| aligned[i].as_slice()).collect();
    // C normalizes weights to sum to 1.0 within each group for cpmx_calc_new,
    // then tracks orieff (= raw sum) separately for createcpmxresult blending.
    let w: Vec<f64> = group.iter().map(|&i| weights[i]).collect();
    let sum: f64 = w.iter().sum();
    let wn: Vec<f64> = if sum > 0.0 { w.iter().map(|v| v / sum).collect() } else { vec![1.0; group.len()] };
    let prof = Profile::from_aligned(&seqs, &wn, &scoring.amino_map, scoring.nalphabets);
    (prof, sum)
}

/// Blend two profiles using C's exact createcpmxresult + creategapfreqresult +
/// createogresult + createfgresult logic (MSalignmm.c lines 283-467).
///
/// The ogcp/fgcp blending handles gap positions specially: at block boundaries
/// (gap→non-gap or non-gap→gap), the value is interpolated from the source
/// profile's nongap_freq. Within a gap block, the value is 0.
fn blend_profiles_exact(
    prof1: &Profile,
    prof2: &Profile,
    eff1: f64,
    eff2: f64,
    gaptable1: &[u8],
    gaptable2: &[u8],
    nalphabets: usize,
) -> Profile {
    let alen = gaptable1.len();
    let mut freqs = vec![vec![0.0f64; nalphabets]; alen];
    let mut nongap_freq = vec![0.0f64; alen + 1]; // C uses alen+1
    let mut ogcp = vec![0.0f64; alen];
    let mut fgcp = vec![0.0f64; alen];

    // createcpmxresult: blend frequency matrices
    {
        let mut p = 0usize;
        for j in 0..alen {
            if gaptable1[j] != b'-' {
                if p < prof1.length {
                    for k in 0..nalphabets.min(prof1.freqs[p].len()) {
                        freqs[j][k] += prof1.freqs[p][k] * eff1;
                    }
                }
                p += 1;
            }
        }
    }
    {
        let mut p = 0usize;
        for j in 0..alen {
            if gaptable2[j] != b'-' {
                if p < prof2.length {
                    for k in 0..nalphabets.min(prof2.freqs[p].len()) {
                        freqs[j][k] += prof2.freqs[p][k] * eff2;
                    }
                }
                p += 1;
            }
        }
    }

    // creategapfreqresult: blend nongap frequencies (C uses alen+1 positions)
    {
        let mut p = 0usize;
        for j in 0..=alen {
            if j < alen && gaptable1[j] == b'-' {
                // gap position: skip
            } else {
                if p < prof1.nongap_freq.len() {
                    nongap_freq[j] += prof1.nongap_freq[p] * eff1;
                }
                p += 1;
            }
        }
    }
    {
        let mut p = 0usize;
        for j in 0..alen {
            if gaptable2[j] == b'-' {
                // gap position: skip
            } else {
                if p < prof2.nongap_freq.len() {
                    nongap_freq[j] += prof2.nongap_freq[p] * eff2;
                }
                p += 1;
            }
        }
    }
    nongap_freq[alen] = 1.0; // C: gapfresult[j] = 1.0 at tail

    // createogresult: blend opening gap counts with block-boundary handling
    blend_og_one_side(&mut ogcp, &prof1.ogcp, &prof1.nongap_freq, gaptable1, eff1, prof1.length);
    blend_og_one_side(&mut ogcp, &prof2.ogcp, &prof2.nongap_freq, gaptable2, eff2, prof2.length);

    // createfgresult: blend closing gap counts with block-boundary handling
    blend_fg_one_side(&mut fgcp, &prof1.fgcp, &prof1.nongap_freq, gaptable1, eff1, prof1.length);
    blend_fg_one_side(&mut fgcp, &prof2.fgcp, &prof2.nongap_freq, gaptable2, eff2, prof2.length);

    // Compute gap_freq from nongap_freq
    let gap_freq: Vec<f64> = nongap_freq[..alen].iter().map(|&nf| (1.0 - nf).max(0.0)).collect();
    let nongap_freq_trimmed = nongap_freq[..alen].to_vec();

    Profile {
        freqs,
        gap_freq,
        nongap_freq: nongap_freq_trimmed,
        ogcp,
        fgcp,
        length: alen,
        nalphabets,
    }
}

/// C's createogresult logic for one side (MSalignmm.c lines 354-378).
fn blend_og_one_side(
    result: &mut [f64],
    ori: &[f64],     // raw opening counts
    gf: &[f64],      // nongap_freq
    gaptable: &[u8],
    eff: f64,
    prof_len: usize,
) {
    let alen = result.len();
    let mut p = 0usize;
    for j in 0..alen {
        if gaptable[j] == b'-' {
            if j == 0 {
                result[j] += 1.0 * eff;
            } else if gaptable[j - 1] != b'-' && p > 0 {
                let gf_val = if p - 1 < gf.len() { gf[p - 1] } else { 1.0 };
                result[j] += gf_val * eff;
            }
        } else {
            if j == 0 || (j > 0 && gaptable[j - 1] != b'-') {
                if p < ori.len() {
                    result[j] += ori[p] * eff;
                }
            }
            p += 1;
        }
    }
}

/// C's createfgresult logic for one side (MSalignmm.c lines 419-439).
fn blend_fg_one_side(
    result: &mut [f64],
    ori: &[f64],     // raw closing counts
    gf: &[f64],      // nongap_freq
    gaptable: &[u8],
    eff: f64,
    prof_len: usize,
) {
    let alen = result.len();
    let mut p = 0usize;
    for j in 0..alen {
        if gaptable[j] == b'-' {
            if j == alen - 1 {
                result[j] += eff;
            } else if gaptable[j + 1] != b'-' {
                let gf_val = if p < gf.len() { gf[p] } else { 1.0 };
                result[j] += gf_val * eff;
            }
        } else {
            if j < alen - 1 && gaptable[j + 1] != b'-' {
                if p < ori.len() {
                    result[j] += ori[p] * eff;
                }
            }
            p += 1;
        }
    }
}

fn group_all_gap_columns(group: &[usize], aligned: &[Vec<u8>], width: usize) -> Vec<bool> {
    let mut all_gap = vec![true; width];
    for &idx in group {
        for (col, &ch) in aligned[idx].iter().enumerate() {
            if ch != b'-' && ch != b'.' {
                all_gap[col] = false;
            }
        }
    }
    all_gap
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_tree::{DistanceMatrix, upgma};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};

    fn check_alignment(result: &MultipleAlignment, original: &[Vec<u8>]) {
        let width = result.width();
        assert!(width > 0);
        for (i, seq) in result.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "seq {i} wrong width: {} vs {width}", seq.len());
            let ungapped: Vec<u8> = seq.iter().filter(|&&c| c != b'-').cloned().collect();
            assert_eq!(ungapped, original[i], "seq {i} residues not preserved");
        }
    }

    #[test]
    fn progressive_two_identical() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![b"ACDEFGHIK".to_vec(), b"ACDEFGHIK".to_vec()];
        let names = vec!["s1".into(), "s2".into()];
        let mut dm = DistanceMatrix::new(2);
        dm.set(0, 1, 0.0);
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        assert_eq!(result.sequences[0], result.sequences[1]);
        check_alignment(&result, &seqs);
    }

    #[test]
    fn progressive_three_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIK".to_vec(),
            b"ACDEFHIK".to_vec(),
            b"ACDHIK".to_vec(),
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];
        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.1); dm.set(0, 2, 0.3); dm.set(1, 2, 0.2);
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        check_alignment(&result, &seqs);
    }

    #[test]
    fn progressive_six_sequences_preserves_residues() {
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
        for i in 0..6 { for j in (i+1)..6 { dm.set(i, j, (j-i) as f64 * 0.1); } }
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        check_alignment(&result, &seqs);
    }
}
