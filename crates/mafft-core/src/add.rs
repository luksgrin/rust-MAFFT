/// Add sequences to an existing alignment.
///
/// Ports C `disttbfast.c`'s `--add` (`-K -I N`) flow: concatenate the
/// existing alignment with the new unaligned sequences, build a guide
/// tree over all (existing + new) sequences, compute `mergeoralign[]`
/// to mark each branch as existing-only ('n', skip), or
/// new-touching ('1'/'2'/'w', do the alignment), then run progressive
/// alignment with the existing alignment columns preserved at 'n'
/// branches.
///
/// The existing input sequences must already be aligned (all the same
/// width). The new sequences are treated as raw (gap chars stripped).

use std::cmp::Ordering;

use rayon::prelude::*;

use mafft_align::{build_local_homology_table, profile_align, AlignOp, GapModel, Profile};
use mafft_tree::{
    addonetip_fragment, scoring_matrix_distance, DistanceMatrix, musclesupg, ClusterMethod,
    ktuple_distance, Topology,
};
use mafft_types::{LocalHomologyTable, ScoringContext};

use crate::progressive::{MergeOrAlign, MultipleAlignment};

/// Add new sequences to an existing alignment.
///
/// `existing` is the already-aligned MSA (with gaps).
/// `new_sequences` are the unaligned sequences to add.
/// `new_names` are their names.
/// `scoring` is the scoring context.
/// `use_fft` controls whether FFT acceleration is used for the merge DPs.
///
/// Returns a new MSA containing all sequences (existing + new), aligned.
pub fn add_sequences(
    existing: &MultipleAlignment,
    new_sequences: &[Vec<u8>],
    new_names: &[mafft_types::SeqName],
    scoring: &ScoringContext,
    use_fft: bool,
) -> MultipleAlignment {
    if new_sequences.is_empty() {
        return existing.clone();
    }
    let n_existing = existing.nseq();
    let n_new = new_sequences.len();
    if n_existing == 0 {
        // No anchor; just align the new sequences from scratch.
        return crate::progressive::progressive_align(
            new_sequences,
            new_names,
            &musclesupg(&compute_ktuple_dm(new_sequences), ClusterMethod::default()),
            scoring,
            use_fft,
            None,
        );
    }

    // 1. Strip common-gap columns from the existing alignment. C does
    //    `commongappick(njob-nadd, seq)` at `disttbfast.c:4319` before
    //    the per-step pairalign loop. With our 30-seq existing (width
    //    598), this keeps the columns where at least one existing
    //    sequence has a non-gap residue, dropping the all-gap columns.
    let stripped_existing = commongappick(&existing.sequences);

    // 2. Concatenate stripped existing + new (raw).
    let mut all_seqs: Vec<Vec<u8>> = stripped_existing;
    for s in new_sequences {
        // C `gappick0` on every added sequence (`disttbfast.c:4315`,
        // `tbfast.c:3192`): only `-` is a gap. `.` is a legal protein
        // residue (`amino_n['.']` = 23) and must survive into the DP.
        let raw: Vec<u8> = s.iter().filter(|&&c| c != b'-').copied().collect();
        all_seqs.push(raw);
    }
    let mut all_names: Vec<mafft_types::SeqName> = existing.names.clone();
    all_names.extend(new_names.iter().cloned());

    // 3. Guide tree on the all-N+nadd sequences. Distance is computed
    //    on the UNGAPPED form (since stripped existing still has
    //    inter-sequence gaps from the alignment, but the topology
    //    should reflect ungapped similarity).
    let ungapped: Vec<Vec<u8>> = all_seqs
        .iter()
        .map(|s| s.iter().filter(|&&c| c != b'-').copied().collect())
        .collect();
    let dm = compute_ktuple_dm(&ungapped);
    let topo = musclesupg(&dm, ClusterMethod::default());

    // 4. Compute mergeoralign[] tagging each branch.
    let mergeoralign = compute_mergeoralign(&topo, n_existing, n_new);

    // 5. Run progressive alignment with mergeoralign-aware skipping
    //    + new-gap propagation to non-active existing rows.
    crate::progressive::progressive_align_with_mergeoralign_n(
        &all_seqs,
        &all_names,
        &topo,
        &mergeoralign,
        scoring,
        use_fft,
        n_existing,
    )
}

/// Add sequences while preserving the existing alignment's column
/// structure (`--add --keeplength`). After running the standard
/// `add_sequences`, deletes any columns that exist only because of
/// new-sequence insertions, restoring the existing alignment to its
/// original width. Columns where ANY new sequence had a residue but ALL
/// existing sequences had gaps get the new-sequence residue dropped
/// (truncated at that position).
pub fn add_sequences_keeplength(
    existing: &MultipleAlignment,
    new_sequences: &[Vec<u8>],
    new_names: &[mafft_types::SeqName],
    scoring: &ScoringContext,
    use_fft: bool,
) -> MultipleAlignment {
    let (msa, _) = add_sequences_keeplength_with_map(
        existing, new_sequences, new_names, scoring, use_fft,
    );
    msa
}

/// Like [`add_sequences_keeplength`] but also returns the per-added-
/// sequence list of dropped insertion runs (`(start_in_addbk_0based,
/// run_length)` pairs). Used by `--mapout` / `--compactmapout`.
///
/// Each entry corresponds to a maximal run of residues in the
/// ORIGINAL (pre-alignment, gap-free) added sequence that the
/// `--keeplength` column filter dropped — these are the
/// "insertions" the user asked for the map of. Position is
/// 0-indexed in `addbk[i]`.
///
/// C parity: matches the `deletelist` populated by
/// `deletenewinsertions_whole` (`addfunctions.c`) that C feeds to
/// `reconstructdeletemap` / `reconstructdeletemap_compact`.
pub fn add_sequences_keeplength_with_map(
    existing: &MultipleAlignment,
    new_sequences: &[Vec<u8>],
    new_names: &[mafft_types::SeqName],
    scoring: &ScoringContext,
    use_fft: bool,
) -> (MultipleAlignment, Vec<Vec<(usize, usize)>>) {
    if new_sequences.is_empty() {
        return (existing.clone(), Vec::new());
    }
    let n_existing = existing.nseq();
    let nadd = new_sequences.len();
    // C `recordoriginalgaps` (`disttbfast.c:3909-3910`): the existing
    // alignment's all-gap columns, put back after the insertions are
    // deleted (`restoreoriginalgaps`, `disttbfast.c:4862`).
    let original_gaps = record_original_gaps(&existing.sequences);

    let mut full = add_sequences(existing, new_sequences, new_names, scoring, use_fft);

    let width = full.sequences.first().map(|s| s.len()).unwrap_or(0);
    let mut keep = vec![false; width];
    for col in 0..width {
        let mut any_existing_residue = false;
        for s in full.sequences.iter().take(n_existing) {
            if let Some(&c) = s.get(col) {
                if c != b'-' {
                    any_existing_residue = true;
                    break;
                }
            }
        }
        keep[col] = any_existing_residue;
    }

    // Build per-added-seq deletelist BEFORE column filtering — we
    // need the post-add (pre-filter) sequence to know which residue
    // of `addbk[i]` lives at each dropped column.
    let mut deletelist: Vec<Vec<(usize, usize)>> = Vec::with_capacity(nadd);
    for i in 0..nadd {
        let aligned = &full.sequences[n_existing + i];
        let mut entries: Vec<(usize, usize)> = Vec::new();
        let mut addbk_pos: usize = 0;
        let mut run_start: usize = 0;
        let mut run_len: usize = 0;
        for (col, &c) in aligned.iter().enumerate() {
            if c == b'-' {
                continue;
            }
            if keep[col] {
                if run_len > 0 {
                    entries.push((run_start, run_len));
                    run_len = 0;
                }
            } else {
                if run_len == 0 { run_start = addbk_pos; }
                run_len += 1;
            }
            addbk_pos += 1;
        }
        if run_len > 0 {
            entries.push((run_start, run_len));
        }
        deletelist.push(entries);
    }

    for s in full.sequences.iter_mut() {
        let mut filtered: Vec<u8> = Vec::with_capacity(original_gaps.len());
        for col in 0..s.len() {
            if keep[col] { filtered.push(s[col]); }
        }
        *s = filtered;
    }
    restore_original_gaps(&mut full.sequences, &original_gaps);

    (full, deletelist)
}

/// Compute `mergeoralign[]` for each branch in the guide tree, mirroring
/// C `disttbfast.c:4282-4304` (`--add` non-profile path).
///
/// C's tag uses `includemember(localmem, addmem)` which returns true iff
/// EVERY member of `localmem` is also in `addmem` (i.e., the subtree is
/// *entirely* composed of new sequences). So:
/// - 'n': neither subtree is entirely new (the merge involves at least
///   one existing-or-mixed branch on each side — could be all-existing or
///   a mix of existing + already-merged-new).
/// - '1': LEFT subtree is entirely new, right has at least one
///   non-new (existing or mixed) member.
/// - '2': RIGHT subtree is entirely new, left has at least one
///   non-new member.
/// - 'w': BOTH subtrees are entirely new (a merge of two new-only
///   subclusters).
///
/// For '1' and '2' merges, C strips common gaps from the *non-all-new*
/// (= mixed) side; that's the side that has accumulated gaps from
/// progressive merges.
pub fn compute_mergeoralign(
    topo: &Topology,
    n_existing: usize,
    _n_new: usize,
) -> Vec<MergeOrAlign> {
    let mut tags: Vec<MergeOrAlign> = Vec::with_capacity(topo.steps.len());
    for step in &topo.steps {
        let left_all_new = step.left.iter().all(|&i| i >= n_existing);
        let right_all_new = step.right.iter().all(|&i| i >= n_existing);
        let tag = match (left_all_new, right_all_new) {
            (false, false) => MergeOrAlign::SkipExisting, // 'n'
            (true, false) => MergeOrAlign::NewLeft,        // '1'
            (false, true) => MergeOrAlign::NewRight,       // '2'
            (true, true) => MergeOrAlign::Wide,             // 'w'
        };
        tags.push(tag);
    }
    tags
}

/// `commongappick`: drop columns where ALL `sequences` are gap chars.
/// Mirrors C `addfunctions.c::commongappick`'s filter behavior.
pub fn commongappick(sequences: &[Vec<u8>]) -> Vec<Vec<u8>> {
    if sequences.is_empty() { return Vec::new(); }
    let width = sequences[0].len();
    let mut keep = vec![false; width];
    for col in 0..width {
        for s in sequences {
            let c = s.get(col).copied().unwrap_or(b'-');
            if c != b'-' {
                keep[col] = true;
                break;
            }
        }
    }
    sequences
        .iter()
        .map(|s| {
            let mut out = Vec::with_capacity(width);
            for col in 0..s.len() {
                if keep[col] { out.push(s[col]); }
            }
            out
        })
        .collect()
}

/// Result of [`add_fragments`].
pub(crate) struct FragmentAdd {
    pub msa: MultipleAlignment,
    /// Per fragment, the `(start_in_addbk, run_len)` runs `-Y` dropped.
    pub deletelist: Vec<Vec<(usize, usize)>>,
    /// The `order` file `addsingle` writes for `--reorder`
    /// (`addsingle.c:3601-3615`): the existing tree's root members
    /// (`ordertable`, `:3459`), each followed by the fragments placed
    /// next to it (`follower[]`, `:3581-3598`), in input order.
    pub order: Vec<usize>,
}

/// `sueff` for `addonetip`, which always uses the mixed-linkage formula
/// with `sueff_global` (`mltaln9.c:14464-14465`).
fn cluster_sueff(cluster: ClusterMethod) -> f64 {
    match cluster {
        ClusterMethod::Mix { sueff } => sueff,
        ClusterMethod::Average => 1.0,
        ClusterMethod::Minimum => 0.0,
    }
}

/// `hat2` / `hat2n` are written with `%.3f` and read back by `addsingle`.
fn hat2_round(d: f64) -> f64 {
    (d * 1000.0).round() / 1000.0
}

/// Pair-phase penalties for `dndpre` and `pairlocalalign -Y`.
///
/// The script passes `-f -2.00 -g -0.100 -h 0.100` (`lgop` / `lexp` /
/// `laof`), independent of the progressive gap open used later by
/// `addsingle -f`. `constants()` then scales those the same way L-INS-i
/// does (`ppenalty = (int)(x*1000 - 0.5)`, then `600/1000`).
struct PairPhase {
    matrix_i32: Vec<Vec<i32>>,
    matrix_f64: Vec<Vec<f64>>,
    gap: GapModel,
    score_offset: f64,
}

fn pair_phase(scoring: &ScoringContext) -> PairPhase {
    let nuc = scoring.seq_type.is_nucleotide();
    let (gap_scale, off_scale) = crate::engine::pair_penalty_scales(nuc);
    let cc_int = |x: f64| ((x * 1000.0) - 0.5) as i32;
    let cc_scale = |p: i32, s: f64| ((s * p as f64) + 0.5) as i32;
    let open = cc_scale(cc_int(-2.00), gap_scale);
    let extend = cc_scale(cc_int(-0.100), gap_scale);
    let offset = cc_scale(cc_int(0.100), off_scale);
    let nscored = scoring.nscoredalphabets;
    let mut matrix_i32 = scoring.substitution_matrix.clone();
    let mut matrix_f64 = scoring.consweight_matrix.clone();
    for i in 0..nscored {
        for j in 0..nscored {
            matrix_i32[i][j] -= offset;
            matrix_f64[i][j] -= offset as f64;
        }
    }
    PairPhase {
        matrix_i32,
        matrix_f64,
        gap: GapModel::new(open as f64, extend as f64),
        score_offset: offset as f64 / 600.0,
    }
}

/// `dndpre -M 2`: `maxdist * (1 - naivepairscorefast / min(self))`,
/// clamped `< 0 → 0` and `> 9 → 9.9`, then the hat2 round-trip.
fn dndpre_dm(seqs: &[Vec<u8>], scoring: &ScoringContext, phase: &PairPhase) -> DistanceMatrix {
    let n = seqs.len();
    let mut dm = DistanceMatrix::new(n);
    let penalty = phase.gap.open as i32;
    for i in 0..n {
        for j in (i + 1)..n {
            let mut d = scoring_matrix_distance(
                &seqs[i], &seqs[j], &phase.matrix_i32, &scoring.amino_map, penalty,
            );
            if d < 0.0 {
                d = 0.0;
            }
            if d > 9.0 {
                d = 9.9;
            }
            dm.set(i, j, hat2_round(d));
        }
    }
    dm
}

/// Homology between the existing rows and fragment `iadd` only.
///
/// `readlocalhomtable_one` stores existing × new, and
/// `dontcalcimportance_firstone` (`addsingle.c`, `SMALLMEMORY`) sets
/// `importance = opt * 0.5` instead of running `calcimportance`.
fn fragment_homology(all: &LocalHomologyTable, norg: usize, iadd: usize) -> LocalHomologyTable {
    let src = norg + iadd;
    let mut table = LocalHomologyTable::new(norg + 1);
    for i in 0..norg {
        for region in all.get(i, src) {
            let mut region = region.clone();
            region.importance = region.opt * 0.5;
            region.rimportance = region.importance;
            table.push(i, norg, region);
        }
        for region in all.get(src, i) {
            let mut region = region.clone();
            region.importance = region.opt * 0.5;
            region.rimportance = region.importance;
            table.push(norg, i, region);
        }
    }
    table
}

/// `--addfragments`: C `addsingle` (`scripts/mafft` default branch,
/// `addsingle -Q 100 -W 6 -O -K -I nadd -g -0.01 -F -f -1.53 -h 0.0 -X 0.1`,
/// plus `-Y` when `--keeplength`).
///
/// Each fragment is placed on its own. A 6-mer distance picks a node in
/// the existing UPGMA tree (`addonetip`, with the ungapped-length gate),
/// one profile DP aligns the fragment to that sister clade (`treebase`,
/// new gaps marked `=`), and the independent gap lists are merged
/// (`makegaplistcompact` / `gaplist2alnx` /
/// `insertgapsbyotherfragments_compact`). Blocks where more than one
/// fragment opened a column are realigned (`dorealignment_tree`).
///
/// `scoring.gap.extend` must already be the scaled `-g -0.01` value.
/// This does not run `dvtditr`. `keeplength` is `-Y`.
pub(crate) fn add_fragments(
    existing: &MultipleAlignment,
    new_sequences: &[Vec<u8>],
    new_names: &[mafft_types::SeqName],
    scoring: &ScoringContext,
    use_fft: bool,
    keeplength: bool,
    cluster: ClusterMethod,
) -> FragmentAdd {
    if new_sequences.is_empty() {
        return FragmentAdd {
            msa: existing.clone(),
            deletelist: Vec::new(),
            order: (0..existing.nseq()).collect(),
        };
    }
    let norg = existing.nseq();
    if norg == 0 {
        let msa = crate::progressive::progressive_align(
            new_sequences,
            new_names,
            &musclesupg(&compute_ktuple_dm(new_sequences), cluster),
            scoring,
            use_fft,
            None,
        );
        let n = new_sequences.len();
        return FragmentAdd { msa, deletelist: vec![Vec::new(); n], order: (0..n).collect() };
    }

    let original_gaps = record_original_gaps(&existing.sequences);
    let stripped = commongappick(&existing.sequences);
    let lenfull = stripped.first().map(|s| s.len()).unwrap_or(0);

    let mut ungapped: Vec<Vec<u8>> = stripped
        .iter()
        .map(|s| s.iter().filter(|&&c| c != b'-').copied().collect())
        .collect();
    let new_raw: Vec<Vec<u8>> = new_sequences
        .iter()
        .map(|s| s.iter().filter(|&&c| c != b'-').copied().collect())
        .collect();
    ungapped.extend(new_raw.iter().cloned());
    // Default `--addfragments` (`scripts/mafft` distance=multi, not the
    // large-input 6-mer fallback): existing-existing distances are
    // `dndpre -M 2` on the aligned rows (`-f -2.00 -g -0.100 -h 0.100`),
    // written through hat2's `%.3f`. Existing-to-fragment distances and
    // the homology table are `pairlocalalign -Y` (`L__align11`,
    // scoreoffset 0) on the gap-stripped sequences.
    let phase = pair_phase(scoring);
    let existing_dm = dndpre_dm(&existing.sequences, scoring, &phase);
    let refs: Vec<&[u8]> = ungapped.iter().map(|s| s.as_slice()).collect();
    let (homology, homdist) = build_local_homology_table(
        &refs,
        &phase.matrix_f64,
        &scoring.amino_map,
        &phase.gap,
        phase.score_offset,
    );
    let existing_topo = if norg == 1 {
        Topology::new(1)
    } else {
        musclesupg(&existing_dm, cluster)
    };
    let nogaplen: Vec<usize> = ungapped.iter().take(norg).map(|s| s.len()).collect();
    let alnlen = aln_len_in_node(&stripped, &existing_topo);

    let mut follows = Vec::with_capacity(new_raw.len());
    let mut placed = Vec::with_capacity(new_raw.len());
    let mut deletelist = Vec::with_capacity(new_raw.len());
    let mut kept_frags: Vec<Vec<u8>> = Vec::with_capacity(new_raw.len());

    for (iadd, frag) in new_raw.iter().enumerate() {
        let dists: Vec<f64> = (0..norg)
            .map(|i| {
                if ungapped[i].is_empty() || frag.is_empty() {
                    3.0
                } else {
                    hat2_round(homdist[i][norg + iadd])
                }
            })
            .collect();
        let constraints = fragment_homology(&homology, norg, iadd);
        let added = addonetip_fragment(
            &existing_topo,
            &dists,
            cluster_sueff(cluster),
            &nogaplen,
            frag.len(),
            &alnlen,
        );
        follows.push(added.neighbor);
        let mut topo = added.topology;
        for step in &mut topo.steps {
            if step.left_length < 0.0 {
                step.left_length = 0.0;
            }
            if step.right_length < 0.0 {
                step.right_length = 0.0;
            }
        }
        let mut rows = stripped.clone();
        rows.push(frag.clone());
        let names: Vec<mafft_types::SeqName> = (0..rows.len()).map(|i| i.to_string().into()).collect();
        let tags = compute_mergeoralign(&topo, norg, 1);
        // `addsingle -l 2.7` sets `constraint = 2`: one `A__align` with
        // the homology table, not the FFT merge. `-d` stores that table
        // as existing-row × this fragment only, and
        // `dontcalcimportance_firstone` sets `importance = opt * 0.5`.
        let aligned = crate::progressive::progressive_align_mergeoralign(
            &rows,
            &names,
            &topo,
            &tags,
            scoring,
            false,
            norg,
            b'=',
            false,
            Some(&constraints),
        );
        if keeplength {
            let mut existing_rows = aligned.sequences[..norg].to_vec();
            let mut fragment = aligned.sequences[norg].clone();
            let width = existing_rows.first().map(|s| s.len()).unwrap_or(0);
            fragment.resize(width, b'-');
            let dropped = delete_new_insertions(&mut existing_rows, &mut fragment);
            for c in &mut fragment {
                if *c == b'=' {
                    *c = b'-';
                }
            }
            deletelist.push(dropped);
            kept_frags.push(fragment);
        } else {
            let (gaplist, had) = gaplist_from_private(&aligned.sequences[..norg], lenfull);
            let mut fragment = aligned.sequences.get(norg).cloned().unwrap_or_default();
            for c in &mut fragment {
                if *c == b'=' {
                    *c = b'-';
                }
            }
            placed.push(PlacedFragment { aligned: fragment, gaplist, had_new_gap: had });
        }
    }

    let sequences = if keeplength {
        let mut seqs = stripped;
        seqs.extend(kept_frags);
        restore_original_gaps(&mut seqs, &original_gaps);
        seqs
    } else {
        combine_fragments(stripped, &placed, lenfull, &existing_topo, &follows, scoring)
    };

    let mut names = existing.names.clone();
    names.extend(new_names.iter().cloned());
    let msa = MultipleAlignment {
        sequences,
        names,
        score: 0.0,
        step_trace: Vec::new(),
        guide_tree: None,
        first_pass_sequences: None,
        distance_matrix: None,
    };

    let ordertable: Vec<usize> = match existing_topo.steps.last() {
        Some(root) if norg > 1 => root.left.iter().chain(root.right.iter()).copied().collect(),
        _ => vec![0],
    };
    let mut follower: Vec<Vec<usize>> = vec![Vec::new(); norg];
    for (iadd, &f) in follows.iter().enumerate() {
        follower[f].push(iadd);
    }
    let mut order = Vec::with_capacity(norg + follows.len());
    for &e in &ordertable {
        order.push(e);
        order.extend(follower[e].iter().map(|&a| norg + a));
    }

    FragmentAdd { msa, deletelist, order }
}

struct PlacedFragment {
    aligned: Vec<u8>,
    gaplist: Vec<i32>,
    had_new_gap: bool,
}

#[derive(Clone)]
struct Block {
    nnewres: i32,
    start: i32,
    end: i32,
}

fn gaplist_from_private(existing_rows: &[Vec<u8>], lenfull: usize) -> (Vec<i32>, bool) {
    let rep = existing_rows.iter().find(|row| row.contains(&b'='));
    match rep {
        Some(row) => (makenewgaplist(row, lenfull), true),
        None => (vec![0; lenfull + 2], false),
    }
}

/// `addsingle.c:775`. Counts of `=` before each non-`=` character, plus a
/// trailing count and a `-1` sentinel. `lenfull` is the pre-merge width.
fn makenewgaplist(seq: &[u8], lenfull: usize) -> Vec<i32> {
    let mut l = vec![0i32; lenfull + 2];
    let mut idx = 0usize;
    let mut a = 0usize;
    loop {
        while a < seq.len() && seq[a] == b'=' {
            a += 1;
            if idx < l.len() {
                l[idx] += 1;
            }
        }
        idx += 1;
        if idx < l.len() {
            l[idx] = 0;
        }
        if a >= seq.len() {
            break;
        }
        a += 1;
    }
    if idx < l.len() {
        l[idx] = -1;
    } else if !l.is_empty() {
        let last = l.len() - 1;
        l[last] = -1;
    }
    l
}

/// `addsingle.c:664`.
fn makegaplistcompact(len: usize, posmap: &[i32], gaps: &[i32]) -> Vec<i32> {
    let mut c = vec![0i32; len + 2];
    let mut prep: i32 = -1;
    for i in 0..(len + 2) {
        let pg = posmap[i] - prep - 1;
        if pg > 0 && gaps[i] > 0 {
            c[i] = if pg < gaps[i] { gaps[i] - pg } else { 0 };
        } else {
            c[i] = gaps[i];
        }
        prep = posmap[i];
    }
    c
}

/// `addsingle.c:691`.
fn gaplist2alnx(len: usize, s: &[u8], gaps: &[i32], posmap: &[i32]) -> Vec<u8> {
    let mut a = Vec::new();
    let mut prevp: i32 = -1;
    let mut gi = 0usize;
    let mut pi = 0usize;
    let mut remaining = len;
    while remaining > 0 {
        remaining -= 1;
        let mut gaplen = gaps[gi];
        gi += 1;
        let pval = posmap[pi];
        pi += 1;
        while gaplen > 0 {
            gaplen -= 1;
            a.push(b'-');
        }
        let pos = prevp + 1;
        let mut posl = pval - pos;
        let mut sp = pos;
        while posl > 0 {
            posl -= 1;
            a.push(s[sp as usize]);
            sp += 1;
        }
        a.push(s[pval as usize]);
        prevp = pval;
    }
    let mut gaplen = gaps[gi];
    let pval = posmap[pi];
    while gaplen > 0 {
        gaplen -= 1;
        a.push(b'-');
    }
    let pos = prevp + 1;
    let mut posl = pval - pos;
    let mut sp = pos;
    while posl > 0 {
        posl -= 1;
        a.push(s[sp as usize]);
        sp += 1;
    }
    a
}

/// `addsingle.c:620`. `=` marks a column another fragment already opened.
fn insertgaps_by_others(len: usize, s: &[u8], gaps: &[i32], posmap: &[i32]) -> Vec<u8> {
    let mut a = Vec::new();
    let mut prevp: i32 = -1;
    let mut sp = 0usize;
    for i in 0..len {
        let mut gaplen = gaps[i];
        let pval = posmap[i];
        let mut pos = prevp + 1;
        while gaplen > 0 {
            gaplen -= 1;
            pos += 1;
            a.push(s[sp]);
            sp += 1;
        }
        while pos < pval {
            pos += 1;
            a.push(b'=');
        }
        a.push(s[sp]);
        sp += 1;
        prevp = pval;
    }
    let mut gaplen = gaps[len];
    let pval = posmap[len];
    let mut pos = prevp + 1;
    while gaplen > 0 {
        gaplen -= 1;
        pos += 1;
        a.push(s[sp]);
        sp += 1;
    }
    while pos < pval {
        pos += 1;
        a.push(b'=');
    }
    a
}

/// `addsingle.c:575`. Updates `posmap[0..=len]`.
fn adjustposmap(len: usize, posmap: &mut [i32], gaps: &[i32]) {
    let mut newpos = vec![0i32; len + 2];
    let mut zure = 0i32;
    for i in 0..len {
        zure += gaps[i];
        newpos[i] = posmap[i] + zure;
    }
    zure += gaps[len];
    newpos[len] = posmap[len] + zure;
    for i in 0..=len {
        posmap[i] = newpos[i];
    }
}

/// `addsingle.c:130`.
fn countnewres(len: usize, blocks: &mut [Block], posmap: &[i32], gaps: &[i32]) {
    let mut regstart: i32 = 0;
    for i in 0..=len {
        if blocks[i].nnewres != 0 || gaps[i] != 0 {
            blocks[i].start = regstart;
            blocks[i].end = posmap[i] - 1;
        }
        if gaps[i] != 0 {
            blocks[i].nnewres += 1;
        }
        regstart = posmap[i] + 1;
    }
}

fn combine_fragments(
    stripped: Vec<Vec<u8>>,
    placed: &[PlacedFragment],
    lenfull: usize,
    existing_topo: &Topology,
    follows: &[usize],
    scoring: &ScoringContext,
) -> Vec<Vec<u8>> {
    let norg = stripped.len();
    let mut posmap: Vec<i32> = (0..=lenfull as i32).collect();
    posmap.push(0);
    let mut blocks = vec![Block { nnewres: 0, start: 0, end: 0 }; lenfull + 1];
    let mut seqs = stripped;
    for (iadd, frag) in placed.iter().enumerate() {
        let compact = makegaplistcompact(lenfull, &posmap, &frag.gaplist);
        if iadd == 0 || frag.had_new_gap {
            for row in &mut seqs {
                *row = gaplist2alnx(lenfull, row, &compact, &posmap);
            }
        }
        let mut added = insertgaps_by_others(lenfull, &frag.aligned, &frag.gaplist, &posmap);
        for c in &mut added {
            if *c == b'=' {
                *c = b'-';
            }
        }
        seqs.push(added);
        adjustposmap(lenfull, &mut posmap, &compact);
        countnewres(lenfull, &mut blocks, &posmap, &frag.gaplist);
    }
    for i in 0..=lenfull {
        if blocks[i].nnewres > 1 {
            let zure = dorealignment_tree(
                &mut seqs, blocks[i].start, blocks[i].end, norg, existing_topo, follows, scoring,
            );
            for j in (i + 1)..=lenfull {
                if blocks[j].nnewres != 0 {
                    blocks[j].start -= zure;
                    blocks[j].end -= zure;
                }
            }
        }
    }
    seqs
}

struct FragPiece {
    seq_index: usize,
    data: Vec<u8>,
}

/// `addsingle.c:175` `dorealignment_tree`. Returns `zure` (old block
/// length minus the realigned length) so later blocks shift.
fn dorealignment_tree(
    seqs: &mut [Vec<u8>],
    start: i32,
    end: i32,
    norg: usize,
    existing_topo: &Topology,
    follows: &[usize],
    scoring: &ScoringContext,
) -> i32 {
    if start < 0 || end < start {
        return 0;
    }
    let width = seqs.first().map(|s| s.len()).unwrap_or(0);
    let start_u = start as usize;
    let end_u = end as usize;
    if end_u >= width {
        return 0;
    }
    let old_len = (end - start + 1) as i32;
    let mut pieces: Vec<FragPiece> = Vec::new();
    for i in 0..seqs.len() {
        let slice = &seqs[i][start_u..=end_u];
        let residues: Vec<u8> = slice.iter().copied().filter(|&c| c != b'-').collect();
        if !residues.is_empty() && i >= norg {
            pieces.push(FragPiece { seq_index: i, data: residues });
        }
    }
    if pieces.is_empty() {
        return 0;
    }

    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut g2n: Vec<usize> = Vec::new();
    for i in 0..pieces.len() {
        let nearest = follows[pieces[i].seq_index - norg];
        let g = groups.iter().position(|members| {
            let n0 = pieces[members[0]].seq_index - norg;
            follows[n0] == nearest
        });
        let g = match g {
            Some(g) => g,
            None => {
                g2n.push(nearest);
                groups.push(Vec::new());
                groups.len() - 1
            }
        };
        groups[g].push(i);
    }
    for group in &mut groups {
        mafft_tree::bsd_qsort::bsd_qsort(group, |a, b| {
            let la = pieces[*a].data.len();
            let lb = pieces[*b].data.len();
            if la > lb {
                Ordering::Less
            } else if la < lb {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });
        for j in 1..group.len() {
            let base: Vec<usize> = group[..j].to_vec();
            let one = [group[j]];
            profile_align_pieces(&mut pieces, &one, &base, scoring);
        }
    }

    let mut pickhistory: Vec<Option<usize>> = vec![None; groups.len()];
    let mut picks: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    for step in &existing_topo.steps {
        let Some(hit1) = first_group_hit(&step.left, &g2n) else { continue };
        let Some(hit2) = first_group_hit(&step.right, &g2n) else { continue };
        let left = expand_pick(&picks, pickhistory[hit1], hit1);
        let right = expand_pick(&picks, pickhistory[hit2], hit2);
        let nhit = picks.len();
        pickhistory[hit1] = Some(nhit);
        pickhistory[hit2] = Some(nhit);
        picks.push((left, right));
    }
    for (left_g, right_g) in &picks {
        let side_a = flatten_groups(left_g, &groups);
        let side_b = flatten_groups(right_g, &groups);
        profile_align_pieces(&mut pieces, &side_a, &side_b, scoring);
    }

    let newlen = pieces[0].data.len();
    let mut by_seq: Vec<Option<usize>> = vec![None; seqs.len()];
    for (i, piece) in pieces.iter().enumerate() {
        by_seq[piece.seq_index] = Some(i);
    }
    for i in 0..seqs.len() {
        let mut block = if let Some(pi) = by_seq[i] {
            let mut b = pieces[pi].data.clone();
            if b.len() < newlen {
                b.resize(newlen, b'-');
            } else {
                b.truncate(newlen);
            }
            b
        } else {
            vec![b'-'; newlen]
        };
        if block.len() < newlen {
            block.resize(newlen, b'-');
        }
        let tail = seqs[i][end_u + 1..].to_vec();
        seqs[i].truncate(start_u);
        seqs[i].extend_from_slice(&block);
        seqs[i].extend_from_slice(&tail);
    }
    old_len - newlen as i32
}

fn first_group_hit(members: &[usize], g2n: &[usize]) -> Option<usize> {
    for &m in members {
        if let Some(k) = g2n.iter().position(|&g| g == m) {
            return Some(k);
        }
    }
    None
}

fn expand_pick(picks: &[(Vec<usize>, Vec<usize>)], hist: Option<usize>, group: usize) -> Vec<usize> {
    match hist {
        None => vec![group],
        Some(p) => {
            let mut v = picks[p].0.clone();
            v.extend_from_slice(&picks[p].1);
            v
        }
    }
}

fn flatten_groups(group_ids: &[usize], groups: &[Vec<usize>]) -> Vec<usize> {
    let mut out = Vec::new();
    for &g in group_ids {
        out.extend_from_slice(&groups[g]);
    }
    out
}

/// `addfunctions.c:36` `profilealignment2`: uniform weights, `outgap=1`,
/// `newgapstr="-"`, `A__align`. A side of a single character against
/// another single character returns unchanged.
fn profile_align_pieces(
    pieces: &mut [FragPiece],
    side_a: &[usize],
    side_b: &[usize],
    scoring: &ScoringContext,
) {
    if side_a.is_empty() || side_b.is_empty() {
        return;
    }
    if pieces[side_a[0]].data.len() == 1 && pieces[side_b[0]].data.len() == 1 {
        return;
    }
    commongappick_pieces(pieces, side_a);
    commongappick_pieces(pieces, side_b);
    let seqs_a: Vec<Vec<u8>> = side_a.iter().map(|&i| pieces[i].data.clone()).collect();
    let seqs_b: Vec<Vec<u8>> = side_b.iter().map(|&i| pieces[i].data.clone()).collect();
    if seqs_a[0].is_empty() && seqs_b[0].is_empty() {
        return;
    }
    let wa = uniform_weights(&seqs_a);
    let wb = uniform_weights(&seqs_b);
    let refs_a: Vec<&[u8]> = seqs_a.iter().map(|s| s.as_slice()).collect();
    let refs_b: Vec<&[u8]> = seqs_b.iter().map(|s| s.as_slice()).collect();
    let pa = Profile::from_aligned(&refs_a, &wa, &scoring.amino_map, scoring.nalphabets);
    let pb = Profile::from_aligned(&refs_b, &wb, &scoring.amino_map, scoring.nalphabets);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let aln = profile_align(&pa, &pb, &scoring.consweight_matrix, &gap, true, true);
    let (new_a, new_b) = apply_block_ops(&seqs_a, &seqs_b, &aln.operations);
    for (k, &i) in side_a.iter().enumerate() {
        pieces[i].data = new_a[k].clone();
    }
    for (k, &i) in side_b.iter().enumerate() {
        pieces[i].data = new_b[k].clone();
    }
}

fn uniform_weights(seqs: &[Vec<u8>]) -> Vec<f64> {
    let allgap: Vec<bool> = seqs.iter().map(|s| s.iter().all(|&c| c == b'-')).collect();
    let count = allgap.iter().filter(|&&g| !g).count();
    let w = if count == 0 { 0.0 } else { 1.0 / count as f64 };
    allgap.into_iter().map(|g| if g { 0.0 } else { w }).collect()
}

fn apply_block_ops(left: &[Vec<u8>], right: &[Vec<u8>], ops: &[AlignOp]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let w1 = left.first().map(|s| s.len()).unwrap_or(0);
    let w2 = right.first().map(|s| s.len()).unwrap_or(0);
    let mut new_l: Vec<Vec<u8>> = vec![Vec::new(); left.len()];
    let mut new_r: Vec<Vec<u8>> = vec![Vec::new(); right.len()];
    let mut c1 = 0usize;
    let mut c2 = 0usize;
    for op in ops {
        match op {
            AlignOp::Match => {
                for (i, row) in left.iter().enumerate() {
                    new_l[i].push(row.get(c1).copied().unwrap_or(b'-'));
                }
                for (i, row) in right.iter().enumerate() {
                    new_r[i].push(row.get(c2).copied().unwrap_or(b'-'));
                }
                c1 += 1;
                c2 += 1;
            }
            AlignOp::Delete => {
                for (i, row) in left.iter().enumerate() {
                    new_l[i].push(row.get(c1).copied().unwrap_or(b'-'));
                }
                for row in &mut new_r {
                    row.push(b'-');
                }
                c1 += 1;
            }
            AlignOp::Insert => {
                for row in &mut new_l {
                    row.push(b'-');
                }
                for (i, row) in right.iter().enumerate() {
                    new_r[i].push(row.get(c2).copied().unwrap_or(b'-'));
                }
                c2 += 1;
            }
        }
    }
    let _ = (w1, w2);
    (new_l, new_r)
}

fn commongappick_pieces(pieces: &mut [FragPiece], idxs: &[usize]) {
    if idxs.is_empty() {
        return;
    }
    let len = pieces[idxs[0]].data.len();
    if len == 0 {
        return;
    }
    let mut keep = vec![false; len];
    for col in 0..len {
        for &i in idxs {
            if pieces[i].data.get(col).copied().unwrap_or(b'-') != b'-' {
                keep[col] = true;
                break;
            }
        }
    }
    if keep.iter().all(|&k| k) {
        return;
    }
    for &i in idxs {
        let old = std::mem::take(&mut pieces[i].data);
        pieces[i].data = old
            .into_iter()
            .enumerate()
            .filter(|(c, _)| keep[*c])
            .map(|(_, ch)| ch)
            .collect();
    }
}

/// `addfunctions.c:1662` `deletenewinsertions_whole_eq` for one fragment.
/// Columns where every existing row is `-` or `=` are dropped. The
/// returned runs are `(start_in_addbk, length)` of residues that lived
/// in those columns.
fn delete_new_insertions(existing: &mut [Vec<u8>], fragment: &mut Vec<u8>) -> Vec<(usize, usize)> {
    let len = existing.first().map(|s| s.len()).unwrap_or(0);
    if len == 0 {
        return Vec::new();
    }
    let mut eqseq = vec![b'o'; len];
    for i in 0..len {
        let all_gap = existing.iter().all(|s| {
            let c = s.get(i).copied().unwrap_or(b'-');
            c == b'-' || c == b'='
        });
        if all_gap {
            eqseq[i] = b'=';
        }
    }
    let mut entries = Vec::new();
    let mut q = 0usize;
    let mut i = 0usize;
    while i < len {
        let mut tmpc = fragment.get(i).copied().unwrap_or(0);
        let qstart = q;
        while i < len && tmpc != b'-' && tmpc != b'=' && eqseq[i] == b'=' {
            i += 1;
            q += 1;
            tmpc = fragment.get(i).copied().unwrap_or(0);
        }
        if q != qstart {
            entries.push((qstart, q - qstart));
        }
        if tmpc != b'-' && tmpc != b'=' {
            q += 1;
        }
        i += 1;
    }
    let mut p = 0usize;
    for col in 0..len {
        if eqseq[col] != b'=' {
            if p != col {
                for row in existing.iter_mut() {
                    row[p] = row[col];
                }
                if col < fragment.len() {
                    fragment[p] = fragment[col];
                }
            }
            p += 1;
        }
    }
    for row in existing.iter_mut() {
        row.truncate(p);
    }
    fragment.truncate(p);
    entries
}

fn record_original_gaps(seqs: &[Vec<u8>]) -> Vec<u8> {
    let Some(first) = seqs.first() else { return Vec::new() };
    let len = first.len();
    let mut out = vec![b'o'; len];
    for i in 0..len {
        let all_gap = seqs.iter().all(|s| s.get(i).copied().unwrap_or(b'-') == b'-');
        if all_gap {
            out[i] = b'-';
        }
    }
    out
}

fn restore_original_gaps(seqs: &mut [Vec<u8>], original: &[u8]) {
    if original.is_empty() {
        return;
    }
    for seq in seqs.iter_mut() {
        let mut tmp = vec![b'-'; original.len()];
        let mut p = 0usize;
        for (j, &mark) in original.iter().enumerate() {
            if mark == b'-' {
                tmp[j] = b'-';
            } else if p < seq.len() {
                tmp[j] = seq[p];
                p += 1;
            }
        }
        *seq = tmp;
    }
}

/// Alignment length of each existing-tree node, root first, on a copy
/// (`addsingle.c:1697`). `norg == 1` stores the single row's width at
/// index 0; there is no join step.
fn aln_len_in_node(seqs: &[Vec<u8>], topo: &Topology) -> Vec<usize> {
    if seqs.len() <= 1 {
        return vec![seqs.first().map(|s| s.len()).unwrap_or(0)];
    }
    let mut copy = seqs.to_vec();
    let mut lens = vec![0usize; topo.steps.len()];
    for i in (0..topo.steps.len()).rev() {
        let mut members = topo.steps[i].left.clone();
        members.extend_from_slice(&topo.steps[i].right);
        commongappick_members(&mut copy, &members);
        if let Some(&m) = members.first() {
            lens[i] = copy[m].len();
        }
    }
    lens
}

fn commongappick_members(seqs: &mut [Vec<u8>], members: &[usize]) {
    if members.is_empty() {
        return;
    }
    let len = seqs[members[0]].len();
    let mut keep = vec![false; len];
    for col in 0..len {
        for &m in members {
            if seqs[m].get(col).copied().unwrap_or(b'-') != b'-' {
                keep[col] = true;
                break;
            }
        }
    }
    if keep.iter().all(|&k| k) {
        return;
    }
    for &m in members {
        let old = std::mem::take(&mut seqs[m]);
        seqs[m] = old
            .into_iter()
            .enumerate()
            .filter(|(c, _)| *c < keep.len() && keep[*c])
            .map(|(_, ch)| ch)
            .collect();
    }
}

fn compute_ktuple_dm(sequences: &[Vec<u8>]) -> DistanceMatrix {
    let nseq = sequences.len();
    let pairs: Vec<(usize, usize, f64)> = (0..nseq)
        .into_par_iter()
        .flat_map(|i| {
            let seqs = sequences;
            ((i + 1)..nseq).into_par_iter().map(move |j| {
                let d = ktuple_distance(&seqs[i], &seqs[j], 6);
                (i, j, d)
            })
        })
        .collect();

    let mut dm = DistanceMatrix::new(nseq);
    for (i, j, d) in pairs {
        dm.set(i, j, d);
    }
    dm
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_scoring::build_context;
    use mafft_tree::JoinStep;
    use mafft_types::{ScoringModel, SeqType};

    fn make_existing_alignment() -> MultipleAlignment {
        MultipleAlignment {
            sequences: vec![
                b"ACDEFGHIK".to_vec(),
                b"ACDEF-HIK".to_vec(),
                b"ACD---HIK".to_vec(),
            ],
            names: vec!["s1".into(), "s2".into(), "s3".into()],
            score: 0.0,
            step_trace: Vec::new(), guide_tree: None, first_pass_sequences: None, distance_matrix: None,
        }
    }

    #[test]
    fn add_single_sequence() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let new_seqs = vec![b"ACDEFGHIKLM".to_vec()];
        let new_names = vec!["new1".into()];
        let result = add_sequences(&existing, &new_seqs, &new_names, &scoring, false);
        assert_eq!(result.nseq(), 4);
        let w = result.sequences[0].len();
        for s in &result.sequences { assert_eq!(s.len(), w); }
    }

    #[test]
    fn add_no_sequences_returns_existing() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let result = add_sequences(&existing, &[], &[], &scoring, false);
        assert_eq!(result.nseq(), 3);
    }

    #[test]
    fn keeplength_preserves_existing_width() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let target_width = existing.sequences[0].len();
        let new_seqs = vec![b"ACDEFGHIK".to_vec()];
        let new_names = vec!["new1".into()];
        let result = add_sequences_keeplength(&existing, &new_seqs, &new_names, &scoring, false);
        assert_eq!(result.nseq(), 4);
        for s in &result.sequences {
            assert_eq!(s.len(), target_width);
        }
    }

    #[test]
    fn mergeoralign_three_existing_one_new() {
        // Topology: 4 leaves, 3 join steps. Existing = [0, 1, 2], new = [3].
        let mut topo = Topology::new(4);
        // Step 0: join 0 and 1 (existing-only, 'n').
        topo.steps.push(JoinStep {
            left: vec![0], right: vec![1], left_length: 0.0, right_length: 0.0,
        });
        // Step 1: join (0,1) with 2 (existing-only, 'n').
        topo.steps.push(JoinStep {
            left: vec![0, 1], right: vec![2], left_length: 0.0, right_length: 0.0,
        });
        // Step 2: join (0,1,2) with 3 (right is all-new, '2').
        topo.steps.push(JoinStep {
            left: vec![0, 1, 2], right: vec![3], left_length: 0.0, right_length: 0.0,
        });

        let tags = compute_mergeoralign(&topo, 3, 1);
        assert_eq!(tags.len(), 3);
        assert!(matches!(tags[0], MergeOrAlign::SkipExisting));
        assert!(matches!(tags[1], MergeOrAlign::SkipExisting));
        assert!(matches!(tags[2], MergeOrAlign::NewRight));
    }

    #[test]
    fn mergeoralign_mixed_subtree_is_n_not_2() {
        // C's `includemember` (mltaln9.c:15053) returns true iff EVERY
        // member of the subtree is in addmem. So a subtree with a MIX
        // of existing and new members is NOT classified as all-new.
        // Topology: 4 leaves; existing = [0,1], new = [2,3].
        let mut topo = Topology::new(4);
        // Step 0: join 2 and 3 (both new, 'w').
        topo.steps.push(JoinStep {
            left: vec![2], right: vec![3], left_length: 0.0, right_length: 0.0,
        });
        // Step 1: join 0 and (2,3) — left is existing, right is all-new ('2').
        topo.steps.push(JoinStep {
            left: vec![0], right: vec![2, 3], left_length: 0.0, right_length: 0.0,
        });
        // Step 2: join (0,2,3) and 1 — left is MIXED, right is existing.
        // Neither side is all-new, so 'n'.
        topo.steps.push(JoinStep {
            left: vec![0, 2, 3], right: vec![1], left_length: 0.0, right_length: 0.0,
        });

        let tags = compute_mergeoralign(&topo, 2, 2);
        assert!(matches!(tags[0], MergeOrAlign::Wide));
        assert!(matches!(tags[1], MergeOrAlign::NewRight));
        assert!(matches!(tags[2], MergeOrAlign::SkipExisting));
    }

    #[test]
    fn gaplist_identity_copies_the_row() {
        let s = b"AC-D".to_vec();
        let len = s.len();
        let gaps = vec![0i32; len + 2];
        let mut posmap: Vec<i32> = (0..=len as i32).collect();
        posmap.push(0);
        assert_eq!(gaplist2alnx(len, &s, &gaps, &posmap), s);
    }

    #[test]
    fn gaplist_inserts_leading_dashes() {
        let s = b"ABC".to_vec();
        let mut gaps = vec![0i32; 5];
        gaps[0] = 2;
        let posmap = vec![0, 1, 2, 3, 0];
        assert_eq!(gaplist2alnx(3, &s, &gaps, &posmap), b"--ABC");
    }

    #[test]
    fn makenewgaplist_counts_eq_runs() {
        let seq = b"A=B";
        let l = makenewgaplist(seq, 2);
        assert_eq!(&l[..4], &[0, 1, 0, -1]);
    }
}
