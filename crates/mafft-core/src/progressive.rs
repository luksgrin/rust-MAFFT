/// Progressive alignment following a guide tree.
///
/// Matches C's `treebase()` from `disttbfast.c`: at each merge step,
/// only group1 and group2 sequences are modified. "Other" sequences
/// are left untouched. Profiles are cached after each merge step
/// (`cpmxhist`) to match C's exact float accumulation order.

// §B.5 — `profile_cache` uses `BTreeMap` (not `HashMap`) for deterministic
// iteration order. The cache is currently `.get()` / `.insert()` / `.remove()`
// only, so HashMap's randomized order isn't an active hazard — but any future
// refactor adding `.iter()` / `.values()` / `.keys()` would silently produce
// non-deterministic alignment output. `BTreeMap` makes that future-safe at
// negligible cost (cache size is bounded by `nseq`, keys are small
// `Vec<usize>`).
use std::collections::BTreeMap;
use mafft_align::{profile_align, pairwise_align11, fft_profile_align, Profile, GapModel, AlignOp, FftAlignParams};
use mafft_tree::{Topology, sequence_weights, compute_distfromtip};
use mafft_types::ScoringContext;

/// C `mltaln9.c::dist2offset`: offset = min(0, dist*0.5 - specificityconsideration).
/// `dist` is `2 * distfromtip` so the result is `min(0, distfromtip - sc)`.
pub(crate) fn dist2offset(dist: f64, sc: f64) -> f64 {
    let v = dist * 0.5 - sc;
    if v > 0.0 { 0.0 } else { v }
}

/// C `mltaln9.c::makedynamicmtx`. Adds `offset * 600` to every substitution
/// score where `offset = dist2offset(2 * distfromtip, sc)`. Negative for
/// shallow merges (close-related), zero for deep merges. Pulls divergent
/// regions apart at shallow merges → wider final alignment.
///
/// C IMPORTANT: `mltaln9.c:15197-15203` SKIPS cells where amino[i] or
/// amino[j] is '-' (gap_idx in our alphabet). We mirror that here — for
/// protein gap_idx = 24, for DNA gap_idx = 24 (= `b'-'` mapped). Without
/// this skip, profile DP cell scores diverge from C on the `--allowshift`
/// per-step path even when input has no gap characters, because the static
/// `amino_dynamicmtx` in C is char-indexed and the unshifted '-' row/col
/// participates in the boundary handling.
pub(crate) fn make_dynamic_matrix(base: &[Vec<f64>], distfromtip: f64, unalign_level: f64, gap_idx: usize) -> Vec<Vec<f64>> {
    let offset = dist2offset(distfromtip * 2.0, unalign_level);
    if offset == 0.0 {
        return base.iter().map(|r| r.clone()).collect();
    }
    // C `mltaln9.c::makedynamicmtx` computes `out[i][j] = in[i][j] + offset * 600`
    // per cell; clang -O3 with FP_CONTRACT=on fuses this into a single FMA.
    // Pre-computing `delta = offset * 600.0` then `v + delta` is two rounded ops
    // and drifts ~1 ULP per cell. See [[project_allowshift_pairwise_fp]] for the
    // BB12003 bisection. Same shape fix as `constraints.rs` `dyn_matrix` build.
    base.iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter().enumerate()
                .map(|(j, &v)| {
                    if i == gap_idx || j == gap_idx { v } else { offset.mul_add(600.0, v) }
                })
                .collect()
        })
        .collect()
}

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
    /// Final progressive guide tree, kept around so callers can serialize
    /// it for `--treeout` (`mltaln9.c::loadtree` and the various
    /// `fixed_musclesupg_*_treeout` variants). Populated by the engine
    /// when alignment completes; `None` for paths that don't track it
    /// (e.g. tests building an `MultipleAlignment` directly).
    pub guide_tree: Option<mafft_tree::Topology>,
    /// Aligned MSA after the FIRST retree pass, BEFORE the final pass
    /// rebuilds the guide tree and re-aligns. C MAFFT's `--parttree`
    /// runs `splittbfast` twice and feeds CALL 1's output (`pre_1`) to
    /// CALL 2 — both for `--reorder` and `--treeout`. We stash `pre_1`
    /// here so the CLI can replay CALL 2's tree generation on the right
    /// input.
    pub first_pass_sequences: Option<Vec<Vec<u8>>>,
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

/// Per-step branch tag for `--add` progressive alignment, mirroring
/// C `disttbfast.c::mergeoralign[]` (lines 4188-4302):
/// - `SkipExisting` ('n'): both subtrees are existing-only — the
///   alignment is already in place; do nothing.
/// - `NewLeft` ('1'): only the LEFT subtree contains a "new" sequence.
/// - `NewRight` ('2'): only the RIGHT subtree contains a "new" sequence.
/// - `Wide` ('w'): both subtrees contain new sequences — full DP merge.
#[derive(Debug, Clone, Copy)]
pub enum MergeOrAlign {
    SkipExisting,
    NewLeft,
    NewRight,
    Wide,
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
        sequences, names, topology, scoring, use_fft, shift_penalty, None, false,
    )
}

/// Run progressive alignment with all per-sequence weights set to 1.0.
/// When normalized within each cluster, this yields uniform weights
/// `1/clus_size` — mirroring C `splittbfast.c::fastconjuction_noweight`'s
/// behavior (used by `--parttree` because `splittbfast.c:6` defines
/// `WEIGHT 0`).
///
/// The standard `progressive_align` derives weights from the guide tree's
/// branch lengths via `weighting::sequence_weights` (matching `disttbfast`
/// / `tbfast`'s `fastconjuction_noname` path).
pub fn progressive_align_unweighted(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
) -> MultipleAlignment {
    let weights = vec![1.0f64; sequences.len()];
    progressive_align_with_weights_override(
        sequences, names, topology, scoring, use_fft, shift_penalty,
        None, false, Some(&weights),
    )
}

/// Progressive alignment for `--add`-style merges, with per-step
/// `mergeoralign[]` tags driving skip/merge decisions. Mirrors C's
/// `disttbfast.c::treebase` mergeoralign-aware loop:
///   - `SkipExisting` branches: do NOTHING — both subtrees are
///     existing-only and already aligned.
///   - All other branches: standard merge_step_cached call, then
///     propagate any new gap columns inserted during the merge to
///     the OTHER existing rows (those not in left/right).
///
/// `sequences[0..n_existing]` should be the (commongappick-stripped)
/// existing alignment, all the same width. `sequences[n_existing..]`
/// are the new raw sequences. The topology is built over all
/// `n_existing + n_new` sequences.
///
/// The gap-propagation step mirrors C's `insertnewgaps_bothorders`
/// (`addfunctions.c:675`) at a coarser grain: instead of using
/// `gaplen`/`gapmap` arrays from `findnewgaps`/`findcommongaps`, we
/// observe the pre-vs-post-merge state of any active existing row to
/// compute new-gap positions and apply the same insertions to the
/// non-active existing rows. Sufficient for byte-identity when the
/// per-step common-gap strip/restore would otherwise be a no-op
/// (= no all-gap columns exist within any active existing subcluster).
pub fn progressive_align_with_mergeoralign(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    mergeoralign: &[MergeOrAlign],
    scoring: &ScoringContext,
    use_fft: bool,
) -> MultipleAlignment {
    progressive_align_with_mergeoralign_n(
        sequences, names, topology, mergeoralign, scoring, use_fft,
        sequences.len(), // default: treat all rows as existing
    )
}

/// Variant that takes `n_existing` explicitly so the caller can
/// distinguish existing rows (whose intra-alignment must be preserved)
/// from new rows (which can be freely re-aligned).
pub fn progressive_align_with_mergeoralign_n(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    mergeoralign: &[MergeOrAlign],
    scoring: &ScoringContext,
    use_fft: bool,
    n_existing: usize,
) -> MultipleAlignment {
    let nseq = sequences.len();
    if nseq == 0 {
        return MultipleAlignment {
            sequences: Vec::new(), names: Vec::new(), score: 0.0, step_trace: Vec::new(), guide_tree: None, first_pass_sequences: None,
        };
    }
    if nseq == 1 {
        return MultipleAlignment {
            sequences: sequences.to_vec(), names: names.to_vec(), score: 0.0, step_trace: Vec::new(), guide_tree: None, first_pass_sequences: None,
        };
    }

    let weights = sequence_weights(topology);
    let mut aligned: Vec<Vec<u8>> = sequences.to_vec();

    let mut last_score = 0.0;
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let mut profile_cache: BTreeMap<Vec<usize>, CachedProfile> = BTreeMap::new();

    // Track which rows are "already aligned" — i.e., have participated
    // in some merge already. Existing rows start aligned (the input
    // alignment). New rows become aligned after their first non-'n'
    // merge. Mirrors C's `alreadyaligned[]` (`disttbfast.c:2257-2258`).
    let mut already_aligned: Vec<bool> = (0..nseq).map(|i| i < n_existing).collect();

    let mut step_trace: Vec<StepTrace> = Vec::with_capacity(topology.steps.len());
    for (step_idx, step) in topology.steps.iter().enumerate() {
        let tag = mergeoralign.get(step_idx).copied().unwrap_or(MergeOrAlign::Wide);
        match tag {
            MergeOrAlign::SkipExisting => {
                let width = aligned[step.left[0]].len().max(aligned[step.right[0]].len());
                step_trace.push(StepTrace {
                    clus1: step.left.len(),
                    clus2: step.right.len(),
                    width,
                    score: 0.0,
                });
            }
            MergeOrAlign::NewRight | MergeOrAlign::NewLeft => {
                // C `disttbfast.c:2745-2756`: for case '2' (NewRight), strip
                // common gaps from group1 (the existing-only side) before the
                // merge, then restore them afterwards. The DP runs between
                // (stripped existing-side) and (full has-new-side). After the
                // merge:
                //   - existing-side and has-new-side are at post_merge_W.
                //   - Restore: re-insert the stripped common-gap columns into
                //     both sides as all-gap columns.
                //   - For OTHER aligned rows (not in either side): expand
                //     from L_pre to L_pre + N1 by inserting gap chars at the
                //     positions corresponding to new merge gaps in the
                //     existing-side representative.
                //
                // Case '1' (NewLeft) is "nai" (never reached) per
                // `disttbfast.c:2934`, but we handle it symmetrically.
                let (existing_grp, new_grp) = match tag {
                    MergeOrAlign::NewRight => (&step.left[..], &step.right[..]),
                    MergeOrAlign::NewLeft => (&step.right[..], &step.left[..]),
                    _ => unreachable!(),
                };

                let pre_width = aligned[existing_grp[0]].len();

                // findcommongaps: find columns where ALL existing_grp rows
                // are gap. These are the columns to strip.
                let pre_classification: Vec<bool> = (0..pre_width)
                    .map(|col| {
                        existing_grp.iter().all(|&i| {
                            let c = aligned[i].get(col).copied().unwrap_or(b'-');
                            c == b'-' || c == b'.'
                        })
                    })
                    .collect();
                let n_gap_cols = pre_classification.iter().filter(|&&b| b).count();

                // Save the pre-strip representative (for OTHER reconstruction).
                let pre_rep_full = aligned[existing_grp[0]].clone();

                // commongappick(existing_grp): strip the all-gap columns.
                let pre_rep_stripped: Vec<u8>;
                if n_gap_cols > 0 {
                    pre_rep_stripped = pre_rep_full
                        .iter()
                        .enumerate()
                        .filter(|(col, _)| !pre_classification[*col])
                        .map(|(_, &c)| c)
                        .collect();
                    for &i in existing_grp {
                        let stripped: Vec<u8> = aligned[i]
                            .iter()
                            .enumerate()
                            .filter(|(col, _)| !pre_classification[*col])
                            .map(|(_, &c)| c)
                            .collect();
                        aligned[i] = stripped;
                    }
                } else {
                    pre_rep_stripped = pre_rep_full.clone();
                }
                let stripped_width = pre_rep_stripped.len();

                // Run the merge. left/right ordering preserved.
                last_score = merge_step_cached(
                    &step.left,
                    &step.right,
                    &mut aligned,
                    &weights,
                    scoring,
                    &gap,
                    use_fft,
                    &mut profile_cache,
                    None,
                    false,
                    false,
                    false, // c_compat off in --add path
                    true,  // --add: pass 0 only, cache valid
                );

                let post_merge_width = aligned[existing_grp[0]].len();
                let post_rep = aligned[existing_grp[0]].clone();

                // Identify "new merge gap" positions in the post-merge
                // existing-side representative. These are post-merge columns
                // where the merge DP inserted a gap into existing_grp's rows
                // (= positions not corresponding to any pre-strip char).
                let new_merge_gap_set = compute_new_merge_gap_set(&pre_rep_stripped, &post_rep);

                // Build mapping: stripped_idx -> pre-merge anchor positions,
                // and gap_cols_before[s] = the gap_col positions sitting
                // between the (s-1)-th and s-th anchor in pre-merge.
                let anchor_positions: Vec<usize> = (0..pre_width)
                    .filter(|&k| !pre_classification[k])
                    .collect();
                debug_assert_eq!(anchor_positions.len(), stripped_width);
                let mut gap_cols_before: Vec<Vec<usize>> =
                    vec![Vec::new(); stripped_width + 1];
                {
                    let mut s = 0usize;
                    for k in 0..pre_width {
                        if pre_classification[k] {
                            gap_cols_before[s].push(k);
                        } else {
                            s += 1;
                        }
                    }
                }

                // restorecommongaps: for both groups, insert n_gap_cols all-gap
                // columns at the right post-merge positions. For each
                // gap_col at pre-merge position k_pre with stripped_idx_after = s,
                // a gap col is inserted right BEFORE the s-th type A position
                // in post-merge (after any preceding type B's).
                if n_gap_cols > 0 {
                    let inserts_per_strip_idx: Vec<usize> =
                        gap_cols_before.iter().map(|v| v.len()).collect();
                    let active: Vec<usize> =
                        step.left.iter().chain(step.right.iter()).copied().collect();
                    for &i in &active {
                        aligned[i] = restore_common_gaps_to_merged_row(
                            &aligned[i],
                            &new_merge_gap_set,
                            &inserts_per_strip_idx,
                        );
                    }
                }

                // For OTHER already-aligned rows: build the post-restore
                // representation. OTHER preserves its pre-merge content at
                // type A (anchor) and type C (gap_col) positions, and gets
                // gap chars at type B (new merge gap) positions.
                let n1 = post_merge_width - stripped_width;
                if n1 > 0 || n_gap_cols > 0 {
                    let active_set: std::collections::HashSet<usize> =
                        step.left.iter().chain(step.right.iter()).copied().collect();
                    let other_indices: Vec<usize> = (0..nseq)
                        .filter(|i| already_aligned[*i] && !active_set.contains(i))
                        .collect();
                    for i in other_indices {
                        let other_pre = aligned[i].clone();
                        // Defensive: only rebuild if OTHER is at pre_width.
                        // Other widths shouldn't happen if invariants hold.
                        if other_pre.len() == pre_width {
                            aligned[i] = build_other_post_restore_row(
                                &other_pre,
                                &anchor_positions,
                                &gap_cols_before,
                                &new_merge_gap_set,
                                post_merge_width,
                            );
                        }
                    }
                }

                // Mark new-side rows as aligned (existing-side was already).
                for &i in new_grp {
                    already_aligned[i] = true;
                }

                let width = aligned[step.left[0]].len();
                let _ = (pre_width, n_gap_cols, stripped_width, post_merge_width, n1);
                step_trace.push(StepTrace {
                    clus1: step.left.len(),
                    clus2: step.right.len(),
                    width,
                    score: last_score,
                });
            }
            MergeOrAlign::Wide => {
                // Both sides have new sequences. C does no per-step strip
                // for case 'w' (`disttbfast.c:2745-2756` only strips for
                // cases '1' and '2'). Just run the merge.
                last_score = merge_step_cached(
                    &step.left,
                    &step.right,
                    &mut aligned,
                    &weights,
                    scoring,
                    &gap,
                    use_fft,
                    &mut profile_cache,
                    None,
                    false,
                    false,
                    false, // c_compat off in --add path
                    true,  // --add: pass 0 only, cache valid
                );

                for &i in step.left.iter().chain(step.right.iter()) {
                    already_aligned[i] = true;
                }

                let width = aligned[step.left[0]].len().max(aligned[step.right[0]].len());
                step_trace.push(StepTrace {
                    clus1: step.left.len(),
                    clus2: step.right.len(),
                    width,
                    score: last_score,
                });
            }
        }
    }

    let max_width = aligned.iter().map(|s| s.len()).max().unwrap_or(0);
    for seq in &mut aligned {
        seq.resize(max_width, b'-');
    }

    MultipleAlignment {
        sequences: aligned, names: names.to_vec(), score: last_score, step_trace,
        guide_tree: None, first_pass_sequences: None,
    }
}

/// Walk pre-strip and post-merge representatives to identify "new merge gap"
/// post-merge positions. A new merge gap is a post-merge column that doesn't
/// correspond to any pre-strip char (i.e., the DP's "insert" op put a gap in
/// existing_grp at this column).
///
/// Lockstep: walk post; for each post char, if it equals the next pre char,
/// advance both. Otherwise, mark it as a new merge gap. This works because
/// the merge preserves residue order (only inserts gap chars), so pre and
/// post agree on residues with `post_merge_W - L_strip` extra gaps inserted.
fn compute_new_merge_gap_set(
    pre_strip: &[u8],
    post_merge: &[u8],
) -> std::collections::HashSet<usize> {
    let mut gap_set = std::collections::HashSet::new();
    let mut p = 0usize;
    for (q, &c) in post_merge.iter().enumerate() {
        if p < pre_strip.len() && pre_strip[p] == c {
            p += 1;
        } else {
            gap_set.insert(q);
        }
    }
    gap_set
}

/// For a merged-group row at post_merge_W chars, expand to post-restore
/// width by inserting `inserts_per_strip_idx[s]` gap chars right BEFORE
/// the s-th type A (anchor) position, plus trailing
/// `inserts_per_strip_idx[L_strip]` gap chars at the end. Mirrors C's
/// `restorecommongaps` (`addfunctions.c:1453`).
fn restore_common_gaps_to_merged_row(
    row: &[u8],
    new_merge_gap_set: &std::collections::HashSet<usize>,
    inserts_per_strip_idx: &[usize],
) -> Vec<u8> {
    let total_inserts: usize = inserts_per_strip_idx.iter().sum();
    let mut out = Vec::with_capacity(row.len() + total_inserts);
    let mut s = 0usize;
    for q in 0..row.len() {
        if !new_merge_gap_set.contains(&q) {
            for _ in 0..inserts_per_strip_idx[s] {
                out.push(b'-');
            }
            out.push(row[q]);
            s += 1;
        } else {
            out.push(row[q]);
        }
    }
    let l_strip = inserts_per_strip_idx.len() - 1;
    for _ in 0..inserts_per_strip_idx[l_strip] {
        out.push(b'-');
    }
    out
}

/// For an OTHER (already-aligned, not-in-merge) row at pre-merge L_pre,
/// build its post-restore representation. OTHER preserves its char at
/// type A (anchor) and type C (gap_col) positions, and gets gap chars at
/// type B (new merge gap) positions. Mirrors C's `insertnewgaps`
/// (`addfunctions.c:445`) but without the profilealignment refinement.
fn build_other_post_restore_row(
    other_pre: &[u8],
    anchor_positions: &[usize],
    gap_cols_before: &[Vec<usize>],
    new_merge_gap_set: &std::collections::HashSet<usize>,
    post_merge_w: usize,
) -> Vec<u8> {
    let total_size = other_pre.len() + new_merge_gap_set.len();
    let mut out: Vec<u8> = Vec::with_capacity(total_size);
    let mut s = 0usize;
    for q in 0..post_merge_w {
        if new_merge_gap_set.contains(&q) {
            out.push(b'-');
        } else {
            for &k_pre in &gap_cols_before[s] {
                out.push(other_pre.get(k_pre).copied().unwrap_or(b'-'));
            }
            out.push(other_pre.get(anchor_positions[s]).copied().unwrap_or(b'-'));
            s += 1;
        }
    }
    let l_strip = anchor_positions.len();
    for &k_pre in &gap_cols_before[l_strip] {
        out.push(other_pre.get(k_pre).copied().unwrap_or(b'-'));
    }
    out
}

/// Run progressive alignment merges 0..n_steps and return the
/// intermediate `aligned[]` state at that point.
///
/// Used by FFI cross-validation tests that need to reproduce a
/// specific step's input profiles without driving the binary or
/// using env-var-controlled file dumps. Each sequence's length in
/// the returned `Vec<Vec<u8>>` matches whatever cluster width it has
/// at step `n_steps` entry (sequences in different clusters have
/// different widths, mirroring C's progressive merge state).
///
/// Passing `n_steps == topology.steps.len()` runs all merges and
/// returns the final padded alignment.
pub fn progressive_align_partial(
    sequences: &[Vec<u8>],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
    n_steps: usize,
) -> Vec<Vec<u8>> {
    let nseq = sequences.len();
    if nseq <= 1 {
        return sequences.to_vec();
    }

    let weights = sequence_weights(topology);
    let mut aligned: Vec<Vec<u8>> = sequences.to_vec();

    let mut gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    if let Some(shift) = shift_penalty {
        gap = gap.with_shift(shift);
    }

    let mut profile_cache: BTreeMap<Vec<usize>, CachedProfile> = BTreeMap::new();

    let limit = n_steps.min(topology.steps.len());
    for step in topology.steps.iter().take(limit) {
        merge_step_cached(
            &step.left, &step.right, &mut aligned, &weights, scoring, &gap, use_fft,
            &mut profile_cache, None, false, false,
            false, // c_compat off in partial replay (test diagnostics)
            true,  // partial replay mirrors pass 0
        );
    }
    aligned
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
    penalize_term_gaps: bool,
) -> MultipleAlignment {
    progressive_align_with_weights_override(
        sequences, names, topology, scoring, use_fft, shift_penalty,
        constraints, penalize_term_gaps, None,
    )
}

/// Like `progressive_align_with_constraints` but allows overriding the
/// per-sequence weights. When `weights_override` is `Some(w)`, each
/// `w[i]` is used directly (still normalized within each cluster at
/// merge time). When `None`, the weights come from
/// `sequence_weights(topology)` (the tree-derived
/// `weightFromABranch`-based defaults).
pub fn progressive_align_with_weights_override(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
    penalize_term_gaps: bool,
    weights_override: Option<&[f64]>,
) -> MultipleAlignment {
    progressive_align_full(
        sequences, names, topology, scoring, use_fft, shift_penalty,
        constraints, penalize_term_gaps, weights_override, 0.0, false, false,
    )
}

/// Like `progressive_align_with_weights_override` but with `unalign_level`
/// (C's `specificityconsideration`). When `unalign_level > 0`, builds a
/// per-step dynamic substitution matrix that scales by `(distfromtip -
/// unalign_level) * 600` (clamped at 0). Mirrors `disttbfast.c:2304-2307` +
/// `mltaln9.c::makedynamicmtx`. `--allowshift` sets this to 0.8.
pub fn progressive_align_full(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
    penalize_term_gaps: bool,
    weights_override: Option<&[f64]>,
    unalign_level: f64,
    legacy_gap_cost: bool,
    memsave_dp: bool,
) -> MultipleAlignment {
    progressive_align_full_c_compat(
        sequences, names, topology, scoring, use_fft, shift_penalty,
        constraints, penalize_term_gaps, weights_override, unalign_level,
        legacy_gap_cost, memsave_dp, false,
    )
}

/// Like `progressive_align_full` but with an explicit `c_compat` flag
/// that enables C MAFFT's static-TLS cpmx memoization (see
/// `mafft_align::profile::CPMX_MEMO`). Default callers should use
/// `progressive_align_full` which passes `c_compat=false`.
pub fn progressive_align_full_c_compat(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
    penalize_term_gaps: bool,
    weights_override: Option<&[f64]>,
    unalign_level: f64,
    legacy_gap_cost: bool,
    memsave_dp: bool,
    c_compat: bool,
) -> MultipleAlignment {
    progressive_align_full_c_compat_ex(
        sequences, names, topology, scoring, use_fft, shift_penalty,
        constraints, penalize_term_gaps, weights_override, unalign_level,
        legacy_gap_cost, memsave_dp, c_compat, true,
    )
}

/// Like `progressive_align_full_c_compat` but with an explicit `use_cache`
/// flag. When `use_cache` is false, every merge rebuilds its child profiles
/// fresh from the current `aligned[]` state via `cpmx_calc_new`, mirroring
/// C MAFFT's `dooneiteration` behavior (`disttbfast.c:2288-2289` —
/// `cpmxchild0/1 = NULL`). When true (default), child profiles are blended
/// from cached parent profiles via `blend_profiles_exact`, mirroring C's
/// `createcpmxresult` in `treebase` (pass 0). The blend matches C's
/// `createcpmxresult` exactly (including its "tsukawanai" comment that
/// excludes the eff*1.0 gap contribution at gap-insertion positions —
/// see `Salignmm.c:622-626`); using it in pass-1+ would diverge from C,
/// which forces fresh `cpmx_calc_new` there.
pub fn progressive_align_full_c_compat_ex(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
    penalize_term_gaps: bool,
    weights_override: Option<&[f64]>,
    unalign_level: f64,
    legacy_gap_cost: bool,
    memsave_dp: bool,
    c_compat: bool,
    use_cache: bool,
) -> MultipleAlignment {
    let nseq = sequences.len();
    if nseq == 0 {
        return MultipleAlignment {
            sequences: Vec::new(), names: Vec::new(), score: 0.0, step_trace: Vec::new(), guide_tree: None, first_pass_sequences: None,
        };
    }
    if nseq == 1 {
        return MultipleAlignment {
            sequences: sequences.to_vec(), names: names.to_vec(), score: 0.0, step_trace: Vec::new(), guide_tree: None, first_pass_sequences: None,
        };
    }

    let weights = match weights_override {
        Some(w) => w.to_vec(),
        None => sequence_weights(topology),
    };
    let mut aligned: Vec<Vec<u8>> = sequences.to_vec();

    let mut last_score = 0.0;
    let mut gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64)
        .with_legacy_gap_cost(legacy_gap_cost);
    if let Some(shift) = shift_penalty {
        gap = gap.with_shift(shift);
    }

    // Profile cache: maps a set of sequence indices (sorted) to its cached profile.
    // After each merge, the merged profile is stored so the next merge can reuse it.
    let mut profile_cache: BTreeMap<Vec<usize>, CachedProfile> = BTreeMap::new();

    // `--c-compat`: reset per-thread cpmx memo at start of this pass
    // (mirrors C `Salignmm.c:1365-1366` which sets previousfirstlen=-1,
    // previousicyc=-1 on buffer resize). For multi-pass progressive
    // (retree>1), each pass starts with a clean memo so cross-pass
    // state doesn't leak.
    if c_compat {
        mafft_align::reset_cpmx_memo();
    }

    // Per-step dynamic-matrix offset. `--allowshift`/`--unalignlevel`
    // triggers `unalign_level > 0`. C builds a fresh `dynamicmtx` per
    // step from the tree node height (`disttbfast.c:2304-2307`).
    let distfromtip: Vec<f64> = if unalign_level > 0.0 {
        compute_distfromtip(topology)
    } else {
        Vec::new()
    };
    // Per-step scoring contexts (only when unalign_level > 0). Each
    // entry differs from `scoring` only in `consweight_matrix` (the
    // f64 version of the substitution matrix used by all DP routines).
    let gap_idx = scoring.amino_map[b'-' as usize] as usize;
    let dyn_scoring: Vec<ScoringContext> = if unalign_level > 0.0 {
        distfromtip
            .iter()
            .map(|&dft| {
                let mut s = scoring.clone();
                s.consweight_matrix =
                    make_dynamic_matrix(&scoring.consweight_matrix, dft, unalign_level, gap_idx);
                s
            })
            .collect()
    } else {
        Vec::new()
    };

    let mut step_trace: Vec<StepTrace> = Vec::with_capacity(topology.steps.len());
    for (step_idx, step) in topology.steps.iter().enumerate() {
        let step_scoring: &ScoringContext = if unalign_level > 0.0 {
            &dyn_scoring[step_idx]
        } else {
            scoring
        };
        last_score = merge_step_cached(
            &step.left, &step.right, &mut aligned, &weights, step_scoring, &gap, use_fft,
            &mut profile_cache, constraints, penalize_term_gaps, memsave_dp, c_compat,
            use_cache,
        );

        let width = aligned[step.left[0]].len().max(aligned[step.right[0]].len());
        step_trace.push(StepTrace {
            clus1: step.left.len(),
            clus2: step.right.len(),
            width,
            score: last_score,
        });
        if std::env::var("MAFFT_DEBUG_STEPS").is_ok() {
            eprintln!("RDBG {} {} {} {} {:.4}",
                step_idx, step.left.len(), step.right.len(), width, last_score);
        }
        if let Ok(f) = std::env::var("RS_PROGRESSIVE_TRACE") {
            use std::io::Write;
            if let Ok(mut fp) = std::fs::OpenOptions::new().create(true).append(true).open(&f) {
                let m1 = step.left[0];
                let m2 = step.right[0];
                let mut h: u64 = 5381;
                for &i in &step.left {
                    for &c in &aligned[i] { h = h.wrapping_mul(33).wrapping_add(c as u64); }
                }
                for &i in &step.right {
                    for &c in &aligned[i] { h = h.wrapping_mul(33).wrapping_add(c as u64); }
                }
                let _ = writeln!(fp, "step={} m1={} m2={} clus1={} clus2={} width={} pscore={:.6} hash={:x}",
                    step_idx, m1, m2, step.left.len(), step.right.len(), width, last_score, h);
            }
        }
        if std::env::var("RDBG_PT_STEPS").is_ok() {
            eprintln!("RDBG_PT step={} clus1={} clus2={} width={} mem1={:?} mem2={:?}",
                step_idx, step.left.len(), step.right.len(), width, step.left, step.right);
        }
        // BB30013 cpmxhist diagnostic: dump the cached profile for this
        // step's output cluster after the merge writes to cache. Matches
        // C's `disttbfast.c::treebase` cpmxhist dump at the same point.
        if let Ok(prefix) = std::env::var("MAFFT_DUMP_CPMX_PREFIX") {
            let mut key = step.left.clone();
            key.extend_from_slice(&step.right);
            key.sort_unstable();
            if let Some(cached) = profile_cache.get(&key) {
                let fname = format!("{}_step_{}.txt", prefix, step_idx);
                if let Ok(mut f) = std::fs::File::create(&fname) {
                    use std::io::Write;
                    let prof = &cached.profile;
                    let cw = prof.length;
                    let na = prof.nalphabets;
                    writeln!(f, "step={} width={} nalphabets={} clus1={} clus2={} score={:.6}",
                        step_idx, cw, na, step.left.len(), step.right.len(), last_score).unwrap();
                    for k in 0..na {
                        write!(f, "F[{}]:", k).unwrap();
                        for j in 0..cw {
                            write!(f, " {:.18e}", prof.freqs[j][k]).unwrap();
                        }
                        writeln!(f).unwrap();
                    }
                    // C's `gapfreq*pt` stores `nongap_freq` (= 1.0 - gap_freq);
                    // see `Salignmm.c:1495,1519` for the post-gapcountf flip.
                    // Rust caches `nongap_freq` of length `cw` and sets
                    // `nongap_freq[cw] = 1.0` implicitly in the DP. To match
                    // C's cpmxhist[nalphabets] which has length cw+1, we
                    // emit cw nongap_freq values + the implied 1.0 terminator.
                    write!(f, "G:").unwrap();
                    for j in 0..cw {
                        write!(f, " {:.18e}", prof.nongap_freq[j]).unwrap();
                    }
                    write!(f, " {:.18e}", 1.0).unwrap();
                    writeln!(f).unwrap();
                    write!(f, "O:").unwrap();
                    for j in 0..cw {
                        write!(f, " {:.18e}", prof.ogcp[j]).unwrap();
                    }
                    writeln!(f).unwrap();
                    write!(f, "N:").unwrap();
                    for j in 0..cw {
                        write!(f, " {:.18e}", prof.fgcp[j]).unwrap();
                    }
                    writeln!(f).unwrap();
                }
            }
        }
    }

    let max_width = aligned.iter().map(|s| s.len()).max().unwrap_or(0);
    for seq in &mut aligned {
        seq.resize(max_width, b'-');
    }

    MultipleAlignment {
        sequences: aligned, names: names.to_vec(), score: last_score, step_trace,
        guide_tree: None, first_pass_sequences: None,
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
    cache: &mut BTreeMap<Vec<usize>, CachedProfile>,
    constraints: Option<&mafft_types::LocalHomologyTable>,
    penalize_term_gaps: bool,
    memsave_dp: bool,
    c_compat: bool,
    use_cache: bool,
) -> f64 {
    let width1 = aligned[group1[0]].len();
    let width2 = aligned[group2[0]].len();

    // Look up cached profiles or build from sequences
    let key1 = sorted_key(group1);
    let key2 = sorted_key(group2);

    // C-compat path: when enabled, build prof1 via `from_aligned_with_memo`
    // so the thread-local CPMX_MEMO incrementally updates when conditions
    // match C's `reuseprofiles` (Salignmm.c:1446-1450). For C, only the
    // FIRST cluster (cluster1) participates in the memo — cluster2 is
    // always rebuilt from scratch (`cpmx_calc_new(seq2, ...)` at
    // Salignmm.c:1555). Mirror that asymmetry here.
    // C MAFFT only uses the cpmxhist cache in pass 0 (`treebase`); pass 1+
    // (`dooneiteration`) sets `cpmxchild0/1 = NULL` and falls through to
    // `cpmx_calc_new` (Salignmm.c:1473-1505 fallback path,
    // disttbfast.c:2288-2289). The cached blend (`createcpmxresult`) omits
    // the eff*1.0 gap-insertion contribution to `cpmx[24][j]` ("tsukawanai"
    // comment at Salignmm.c:624), so reusing it across passes diverges from
    // a fresh build — this surfaces as the BB20018 / BB40046 step-50
    // pass-1 divergences. `use_cache=false` here forces fresh
    // `from_aligned` rebuilds in refinement passes to match C.
    let try_cache = use_cache;
    let (prof1, eff1) = if try_cache && cache.get(&key1).is_some() {
        let cached = cache.get(&key1).unwrap();
        (cached.profile.clone(), cached.eff)
    } else if c_compat {
        build_profile_with_memo(group1, aligned, weights, scoring)
    } else {
        let (prof, eff) = build_profile_from_seqs(group1, aligned, weights, scoring);
        (prof, eff)
    };

    let (prof2, eff2) = if try_cache && cache.get(&key2).is_some() {
        let cached = cache.get(&key2).unwrap();
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
        // G__align11 path: flat gap penalty, character-level scoring.
        // Only used when FFT is disabled and no constraints. With constraints
        // (L-INS-i / E-INS-i / G-INS-i), single-vs-single merges still need
        // the impmtx contribution per cell — fall through to the constrained
        // profile DP below.
        //
        // `head_gap` / `tail_gap` must mirror C's `outgap` (= `penalize_term_gaps`).
        // For `--parttree --nofft`, the script omits `-O` so `outgap=1`
        // (`scripts/mafft:2655`, `splittbfast.c:560`), meaning every per-pair
        // merge — including 1-vs-1 — must penalize terminal gaps. Hardcoding
        // `false` here caused the `--parttree --nofft` 944-line divergence vs C
        // (every 1-vs-1 merge took the wrong head/tail gap path).
        pairwise_align11(
            &aligned[group1[0]], &aligned[group2[0]],
            &scoring.consweight_matrix, &scoring.amino_map,
            scoring.gap.open as f64, penalize_term_gaps, penalize_term_gaps,
        )
    } else if use_fft {
        // C uses Falign for ALL steps when use_fft=true (ffttry = nlen > clus,
        // always true). No minimum profile length check.
        //
        // For protein scoring matrices C uses 2-channel polarity+volume FFT
        // via `seq_vec_2` (`Falign.c:342-348`). Build per-internal-index
        // polarity/volume vectors so `find_fft_anchors` can mirror that.
        let property_channels = if scoring.seq_type.is_nucleotide() {
            None
        } else {
            let nscored = scoring.nscoredalphabets;
            let mut polarity_by_idx = vec![0.0f64; nscored];
            let mut volume_by_idx = vec![0.0f64; nscored];
            for ch in 0u16..256 {
                let idx = scoring.amino_map[ch as usize] as usize;
                if idx < nscored {
                    polarity_by_idx[idx] = scoring.polarity[ch as usize];
                    volume_by_idx[idx] = scoring.volume[ch as usize];
                }
            }
            Some((polarity_by_idx, volume_by_idx))
        };
        let fft_params = FftAlignParams {
            num_candidates: 20,
            segment_params: if scoring.seq_type.is_nucleotide() {
                mafft_fft::SegmentParams::dna()
            } else {
                mafft_fft::SegmentParams::protein()
            },
            gap: gap.clone(),
            // C `Falign.c:686-687` sets the per-segment `headgp/tailgp`
            // to the global `outgap` for the first/last segment. So
            // when outgap=1 (term gaps penalized — G-INS-i, --parttree),
            // we set head_gap/tail_gap=true. Threaded via
            // `penalize_term_gaps`.
            head_gap: penalize_term_gaps,
            tail_gap: penalize_term_gaps,
            num_channels: scoring.nscoredalphabets,
            property_channels,
        };
        fft_profile_align(&prof1, &prof2, &scoring.consweight_matrix, &fft_params)
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
        if std::env::var_os("RUST_IMP_DUMP").is_some() {
            let s00 = imp.first().and_then(|r| r.first()).copied().unwrap_or(0.0);
            let s100 = if imp.len()>100 && imp[100].len()>100 { imp[100][100] } else { 0.0 };
            let s300 = if imp.len()>300 && imp[300].len()>300 { imp[300][300] } else { 0.0 };
            eprintln!("[RUST_IMP] g1={:?} g2={:?} lgth1={} lgth2={} eff1={:?} eff2={:?} imp[0,0]={:.4} imp[100,100]={:.4} imp[300,300]={:.4}",
                group1, group2, prof1.length, prof2.length, w1n, w2n, s00, s100, s300);
            for &gi in group1 {
                for &gj in group2 {
                    let regs = table.get(gi, gj);
                    for (idx, r) in regs.iter().enumerate().take(3) {
                        eprintln!("[RUST_IMP] lh[{},{}] e{}: opt={:.6} imp={:.6} overlapaa={} s1={} e1={} s2={} e2={}",
                            gi, gj, idx, r.opt, r.importance, r.overlapaa, r.start1, r.end1, r.start2, r.end2);
                    }
                }
            }
        }
        mafft_align::profile_align_imp(
            &prof1, &prof2, &scoring.consweight_matrix, gap,
            penalize_term_gaps, penalize_term_gaps, Some(&imp),
        )
    } else if memsave_dp {
        // `--memsave`: route through the Hirschberg DP. Mirrors C
        // MAFFT's `MSalignmm` (`tbfast.c:1159-1161` under `alg='M'`).
        // For inputs that fit in memory, the alignment is the same
        // as `profile_align` — only memory layout differs.
        mafft_align::msalignmm(
            &prof1, &prof2, &scoring.consweight_matrix, gap,
            penalize_term_gaps, penalize_term_gaps,
        )
    } else {
        // Non-FFT, no-constraints fallback (`--nofft` path or single-vs-
        // single without constraints). C's `outgap` flows through here
        // via `penalize_term_gaps`: false → outgap=0 (term-gap free,
        // FFT-NS-2/L-INS-i/E-INS-i defaults), true → outgap=1
        // (G-INS-i and `--parttree`).
        profile_align(
            &prof1, &prof2, &scoring.consweight_matrix, gap,
            penalize_term_gaps, penalize_term_gaps,
        )
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

/// `--c-compat` variant of `build_profile_from_seqs` that uses
/// `Profile::from_aligned_with_memo` so the thread-local cpmx memo
/// can fire when conditions match. `firstmem` is `group[0]` (the
/// global leaf index of the cluster's first member, matching C's
/// `localmem[0][0]`); `icyc` is the cluster size; `lgth` is the
/// per-sequence width.
fn build_profile_with_memo(
    group: &[usize],
    aligned: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
) -> (Profile, f64) {
    let seqs: Vec<&[u8]> = group.iter().map(|&i| aligned[i].as_slice()).collect();
    let w: Vec<f64> = group.iter().map(|&i| weights[i]).collect();
    let sum: f64 = w.iter().sum();
    let wn: Vec<f64> = if sum > 0.0 { w.iter().map(|v| v / sum).collect() } else { vec![1.0; group.len()] };
    let firstmem = group[0] as i32;
    let icyc = group.len();
    let lgth = seqs.first().map_or(0, |s| s.len());
    let prof = Profile::from_aligned_with_memo(
        &seqs, &wn, &scoring.amino_map, scoring.nalphabets,
        firstmem, icyc, lgth,
    );
    (prof, sum)
}

/// Blend two profiles using C's exact createcpmxresult + creategapfreqresult +
/// createogresult + createfgresult logic (MSalignmm.c lines 283-467).
///
/// The ogcp/fgcp blending handles gap positions specially: at block boundaries
/// (gap→non-gap or non-gap→gap), the value is interpolated from the source
/// profile's nongap_freq. Within a gap block, the value is 0.
/// Blend two profiles. Exposed `pub` for FFI cross-validation against
/// C's `createcpmxresult + creategapfreqresult + createogresult +
/// createfgresult` (`Salignmm.c:608-823`).
pub fn blend_profiles_exact(
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

    // createcpmxresult: blend frequency matrices.
    // FMA throughout: matches gcc's `-O3` fusion of `a + b*c` so the blended
    // child-profile is bit-identical to C's. Without FMA the per-column
    // weighted frequency accumulates 1-ULP differences that propagate into
    // match_calc_row and flip DP tie-breaks for flat-landscape matrices
    // (TM PAM 200 — §B.2).
    {
        let mut p = 0usize;
        for j in 0..alen {
            if gaptable1[j] != b'-' {
                if p < prof1.length {
                    for k in 0..nalphabets.min(prof1.freqs[p].len()) {
                        freqs[j][k] = prof1.freqs[p][k].mul_add(eff1, freqs[j][k]);
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
                        freqs[j][k] = prof2.freqs[p][k].mul_add(eff2, freqs[j][k]);
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
                    nongap_freq[j] = prof1.nongap_freq[p].mul_add(eff1, nongap_freq[j]);
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
                    nongap_freq[j] = prof2.nongap_freq[p].mul_add(eff2, nongap_freq[j]);
                }
                p += 1;
            }
        }
    }
    nongap_freq[alen] = 1.0; // C: gapfresult[j] = 1.0 at tail

    // createogresult: blend opening gap counts with block-boundary handling
    blend_og_one_side(&mut ogcp, &prof1.ogcp, &prof1.nongap_freq, gaptable1, eff1);
    blend_og_one_side(&mut ogcp, &prof2.ogcp, &prof2.nongap_freq, gaptable2, eff2);

    // createfgresult: blend closing gap counts with block-boundary handling
    blend_fg_one_side(&mut fgcp, &prof1.fgcp, &prof1.nongap_freq, gaptable1, eff1);
    blend_fg_one_side(&mut fgcp, &prof2.fgcp, &prof2.nongap_freq, gaptable2, eff2);

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
) {
    let alen = result.len();
    let mut p = 0usize;
    for j in 0..alen {
        if gaptable[j] == b'-' {
            if j == 0 {
                result[j] += eff;
            } else if gaptable[j - 1] != b'-' && p > 0 {
                let gf_val = if p - 1 < gf.len() { gf[p - 1] } else { 1.0 };
                result[j] = gf_val.mul_add(eff, result[j]);
            }
        } else {
            if j == 0 || (j > 0 && gaptable[j - 1] != b'-') {
                if p < ori.len() {
                    result[j] = ori[p].mul_add(eff, result[j]);
                }
            }
            p += 1;
        }
    }
}

/// C's createfgresult logic for one side (Salignmm.c lines 782-823).
///
/// IMPORTANT: C reads `gaptable1[j+1]` even at j = alen-1, accessing
/// the null terminator past the end of the C string (which is `!= '-'`).
/// So at the LAST non-gap position, C ALWAYS adds `ori[p] * eff` (the
/// closing-count contribution), since the "next" position is treated as
/// non-gap. Our Rust gaptable is a `&[u8]` without a null terminator,
/// so we treat the out-of-bounds j+1 as non-gap to match. The prior
/// `j < alen - 1` short-circuit dropped this contribution at the
/// alignment's last position, causing pass-1 BB20027 cpmx drift when
/// the merged profile was cached and reused at the next merge.
fn blend_fg_one_side(
    result: &mut [f64],
    ori: &[f64],     // raw closing counts
    gf: &[f64],      // nongap_freq
    gaptable: &[u8],
    eff: f64,
) {
    let alen = result.len();
    let mut p = 0usize;
    // C treats out-of-bounds `gaptable[alen]` as non-gap (the `\0` of
    // the null-terminated string, which is != '-').
    let next_is_gap = |j: usize| -> bool {
        if j + 1 < alen { gaptable[j + 1] == b'-' } else { false }
    };
    for j in 0..alen {
        if gaptable[j] == b'-' {
            if j == alen - 1 {
                result[j] += eff;
            } else if !next_is_gap(j) {
                let gf_val = if p < gf.len() { gf[p] } else { 1.0 };
                result[j] = gf_val.mul_add(eff, result[j]);
            }
        } else {
            if !next_is_gap(j) {
                if p < ori.len() {
                    result[j] = ori[p].mul_add(eff, result[j]);
                }
            }
            p += 1;
        }
    }
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
