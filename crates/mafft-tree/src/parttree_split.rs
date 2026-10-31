//! First-level building blocks of `splittbfast.c::splitseq_mq` for
//! `--parttree`, and the public PartTree helpers (topology, `--reorder`
//! order, `--treeout` string) built on the full recursive port in
//! [`crate::parttree_splitseq`].
//!
//! The building blocks here (`build_dfromc`, `assign_to_yukos`,
//! `assemble_topology`, `run_parttree_pipeline_with_scorer`) cover one
//! recursion level: they assume every sequence's yuko group needs no
//! further splitting. That holds only while all pivots survive the
//! redundancy filter, i.e. for fewer than `picksize` (50) sequences
//! (`splittbfast.c:2753-2754` sets `tokyoripara = 0` then). From 50
//! sequences on, yuko groups hold distinct sequences that C aligns in a
//! recursive `splitseq_mq` call before the parent merges them.

use crate::distance::DistanceMatrix;
use crate::musclesupg::{musclesupg, ClusterMethod};
use crate::parttree_pivot::{PartTreePivots, PtSeqKind};
use crate::parttree_dist::{lenfac, common_sextets_p, composition_table, MAX6DIST,
    PLENFACA, PLENFACB, PLENFACC, PLENFACD,
    DLENFACA, DLENFACB, DLENFACC, DLENFACD};
use crate::topology::{JoinStep, Topology};

/// Build `dfromc[i][j]` = parttree distance from the `i`-th surviving
/// yuko to the `j`-th sorted-position sequence (`j` in `[0..nin)`).
/// Mirrors `splittbfast.c:2132-2236`.
pub fn build_dfromc(pivots: &PartTreePivots, kind: PtSeqKind) -> Vec<Vec<f64>> {
    let nyuko = pivots.yukos.len();
    let nin = pivots.scores.len();
    let tsize = match kind { PtSeqKind::Protein => 46656, PtSeqKind::Dna => 4096 };
    let (a, b, c, d) = match kind {
        PtSeqKind::Protein => (PLENFACA, PLENFACB, PLENFACC, PLENFACD),
        PtSeqKind::Dna => (DLENFACA, DLENFACB, DLENFACC, DLENFACD),
    };

    let mut dfromc = vec![vec![0.0f64; nin]; nyuko];

    for (yi, &y_pick) in pivots.yukos.iter().enumerate() {
        let y_sorted_idx = pivots.picks[y_pick]; // sorted-scores index
        let y_score = &pivots.scores[y_sorted_idx];
        let y_table = composition_table(&y_score.points, tsize);

        for j in 0..nin {
            // Self-distance is 0 for the yuko itself.
            if j == y_sorted_idx {
                dfromc[yi][j] = 0.0;
                continue;
            }
            let common = common_sextets_p(&y_table, &pivots.scores[j].points, tsize);
            let bunbo = (y_score.selfscore.min(pivots.scores[j].selfscore)) as f64;
            let raw = if bunbo > 0.0 { 1.0 - common as f64 / bunbo } else { 1.0 };
            let lf = lenfac(y_score.orilen, pivots.scores[j].orilen, a, b, c, d);
            let mut dist = raw * lf;
            if dist > MAX6DIST { dist = MAX6DIST; }
            dfromc[yi][j] = dist;
        }
    }
    dfromc
}

/// Assign every sorted-position sequence (incl. pivots) to its closest
/// yuko. Returns `outs[yi]` = list of original `numinseq` indices
/// assigned to yuko `yi`. Mirrors `splittbfast.c:2244-2302`.
///
/// **Tie-break**: C uses strict `<` (`splittbfast.c:2270-2286`), so the
/// FIRST-encountered yuko (smallest `yi`) wins on ties. Each pivot is
/// closest to itself with distance 0, so it always lands in its own
/// `outs[]` slot.
pub fn assign_to_yukos(
    pivots: &PartTreePivots,
    dfromc: &[Vec<f64>],
) -> Vec<Vec<usize>> {
    let nyuko = pivots.yukos.len();
    let nin = pivots.scores.len();
    let mut outs: Vec<Vec<usize>> = vec![Vec::new(); nyuko];

    for j in 0..nin {
        let mut belongto = 0usize;
        let mut minscore = f64::INFINITY;
        for yi in 0..nyuko {
            if dfromc[yi][j] < minscore {
                minscore = dfromc[yi][j];
                belongto = yi;
            }
        }
        outs[belongto].push(pivots.scores[j].numinseq);
    }
    outs
}

/// Build a yukomtx-as-DistanceMatrix (full square form) for feeding to
/// `musclesupg`, using the `with-diagonal` half-matrix layout that
/// `pivots.yukomtx` uses.
fn yukomtx_to_distance_matrix(pivots: &PartTreePivots) -> DistanceMatrix {
    let nyuko = pivots.yukos.len();
    let mut dm = DistanceMatrix::new(nyuko);
    for i in 0..nyuko {
        for j in (i + 1)..nyuko {
            // pivots.yukomtx[i][j-i] (with-diagonal layout, slot 0 unused).
            dm.set(i, j, pivots.yukomtx[i][j - i]);
        }
    }
    dm
}

/// Assemble the final `Topology` for `--parttree`. Returns one
/// `JoinStep` per UPGMA-on-yukomtx join (`nyuko - 1` steps total).
///
/// Mirrors `splittbfast.c:2444-2519`'s parent pairalign loop: at each
/// UPGMA step, the left/right groups are unions of `outs[yuko_idx]`
/// for the merging yuko-clusters, sorted ascending by original
/// sequence index (matching C's `qsort(mem1, ..., intcompare)` at
/// `splittbfast.c:2477-2478`).
///
/// First recursion level only: a multi-member yuko is emitted as one
/// unaligned group. C instead runs `splitseq_mq` on that group
/// (`splittbfast.c:2377`), which aligns it unless all its members are
/// identical (`nyuko == 1`, `:2105-2123`). Use
/// [`crate::parttree_splitseq::run_splitseq`] for the full algorithm.
pub fn assemble_topology(
    pivots: &PartTreePivots,
    outs: &[Vec<usize>],
) -> Topology {
    let nseq_total: usize = outs.iter().map(|v| v.len()).sum();
    let mut topo = Topology::new(nseq_total);

    // Run UPGMA on the yukomtx and convert each join step into a
    // sequence-level JoinStep. UPGMA's `left` / `right` are sets of
    // yuko-indices; we expand them via `outs[]`.
    let yuko_dm = yukomtx_to_distance_matrix(pivots);
    let yuko_topo = musclesupg(&yuko_dm, ClusterMethod::Mix { sueff: 0.1 });

    for yuko_step in &yuko_topo.steps {
        let mut left_seqs: Vec<usize> = yuko_step.left.iter()
            .flat_map(|&yi| outs[yi].iter().copied())
            .collect();
        let mut right_seqs: Vec<usize> = yuko_step.right.iter()
            .flat_map(|&yi| outs[yi].iter().copied())
            .collect();
        left_seqs.sort();
        right_seqs.sort();
        topo.steps.push(JoinStep {
            left: left_seqs,
            right: right_seqs,
            left_length: yuko_step.left_length,
            right_length: yuko_step.right_length,
        });
    }

    topo
}

/// `--parttree` guide topology: the merge steps of one k-tuple
/// `splittbfast` call, from the recursive `splitseq_mq` port
/// ([`crate::parttree_splitseq`]).
pub fn build_parttree_topology(
    sequences: &[Vec<u8>],
    kind: PtSeqKind,
    picksize: usize,
) -> Topology {
    run_ktuple_with_picksize(sequences, kind, picksize).topology
}

fn run_ktuple_with_picksize(
    sequences: &[Vec<u8>],
    kind: PtSeqKind,
    picksize: usize,
) -> crate::parttree_splitseq::SplitRun {
    use crate::parttree_splitseq::{run_splitseq, KtupleScorer, SplitParams};
    let scorer = KtupleScorer::new(sequences, kind);
    let mut params = SplitParams::for_njob(sequences.len(), false);
    params.picksize = picksize;
    if picksize > sequences.len() {
        params.tokyoripara = 0.0;
    }
    run_splitseq(scorer.entries(), &scorer, params)
}

/// Recursive helper for `c_normalized_yuko_order`: produce the leaf-yuko
/// order for the subtree whose leaves are `leaves`, applying C's
/// smaller-first-element normalization at each internal node (mirroring
/// `mltaln9.c:8184-8197` in `fixed_musclesupg_double_realloc_nobk_halfmtx`).
pub(crate) fn c_normalized_subtree(steps: &[crate::topology::JoinStep], leaves: &[usize]) -> Vec<usize> {
    if leaves.len() == 1 {
        return vec![leaves[0]];
    }
    // Find the join step whose `left ∪ right == leaves`. Search from the
    // back since later steps merge larger sets first (the root is the
    // very last step), and we typically recurse downward.
    let leaf_set: std::collections::BTreeSet<usize> = leaves.iter().copied().collect();
    for step in steps.iter().rev() {
        let combined: std::collections::BTreeSet<usize> =
            step.left.iter().chain(step.right.iter()).copied().collect();
        if combined == leaf_set {
            let l_order = c_normalized_subtree(steps, &step.left);
            let r_order = c_normalized_subtree(steps, &step.right);
            // C: when extending `topol[k][i]`, the two child sub-arrays
            // are concatenated smaller-first-element first
            // (`mltaln9.c:8184-8197`).
            let (first, second) = if !l_order.is_empty()
                && !r_order.is_empty()
                && l_order[0] > r_order[0]
            {
                (r_order, l_order)
            } else {
                (l_order, r_order)
            };
            let mut result = first;
            result.extend(second);
            return result;
        }
    }
    // Fall back to insertion order (shouldn't happen if topology is
    // complete and leaves came from it).
    leaves.to_vec()
}

/// Build the parttree-specific `--treeout` output of one k-tuple
/// `splittbfast` call (`splittbfast.c:1275-1301` leaf emission,
/// `:2532-2553` per-merge concatenation, `:3035` trailing newline).
///
/// Output format differs from the standard musclesupg `--treeout`:
/// - No branch lengths.
/// - Leaves are 1-indexed sequence numbers, NOT sanitized FASTA names.
/// - Multi-member leaf groups render as `(s1,s2,...)` on one line;
///   single-member leaves render as `<n>` sandwiched by `\n`.
/// - No trailing `;` (parttree's tree file omits it).
pub fn compute_parttree_newick(
    sequences: &[Vec<u8>],
    kind: PtSeqKind,
    picksize: usize,
) -> String {
    let nseq = sequences.len();
    if nseq == 0 { return "\n".to_string(); }
    if nseq == 1 { return "\n1\n".to_string(); }
    let mut out = run_ktuple_with_picksize(sequences, kind, picksize).tree;
    out.push('\n');
    out
}

/// The `--reorder` order of one k-tuple `splittbfast` call: the leaf
/// emission order of the recursive `splitseq_mq`
/// (`splittbfast.c:1305-1309`).
pub fn compute_parttree_order(
    sequences: &[Vec<u8>],
    kind: PtSeqKind,
    picksize: usize,
) -> Vec<usize> {
    let nseq = sequences.len();
    if nseq <= 1 {
        return (0..nseq).collect();
    }
    run_ktuple_with_picksize(sequences, kind, picksize).order
}

// =============================================================================
// CALL 2 — second-pass `splittbfast` with `fromaln=1` (`-Z`) scoring.
//
// C MAFFT runs `splittbfast` twice for `--parttree` (`scripts/mafft:2655`
// and `:2681`). The second call uses `-Z` so distances are computed via
// `naivepairscore11` on the already-aligned sequences instead of via
// 6-mer composition. The full recursion lives in `parttree_splitseq`;
// the public helpers below delegate to it.
// =============================================================================

/// Top-level result of a scorer-driven `splitseq_mq` pass: the
/// post-`dcompare` sorted scores, the yuko UPGMA topology, and the
/// per-yuko member lists (`outs[]`) — the first recursion level only.
pub struct PartTreeFromalnResult {
    pub scores: Vec<crate::parttree_pivot::ScoreEntry>,
    pub outs: Vec<Vec<usize>>,
    pub yuko_topo: crate::topology::Topology,
}

/// Run the first level of the splittbfast pivot/yuko/UPGMA pipeline using
/// caller-supplied scoring callbacks. It does not recurse into the yuko
/// groups nor apply the redundancy filter, so it describes C's output only
/// while every pivot survives (fewer than `picksize` sequences); use
/// [`crate::parttree_splitseq::run_splitseq`] for the full algorithm.
///
/// Callbacks operate on the original input indices (0..nseq):
/// - `selfscore`: returns the self-vs-self score for sequence `i`.
/// - `orilen`: returns the ungapped length of sequence `i`.
/// - `pair_score`: returns the pairwise score for `(i, j)`. Implementations
///   compute this however they like (e.g. `naivepairscore11` or
///   `G__align11_noalign`). For `i == j`, returns the selfscore.
/// - `seqs_equal`: returns `true` when `i` and `j` refer to byte-identical
///   sequences (used by the shimon-style dedupe in pivot selection).
pub fn run_parttree_pipeline_with_scorer<S, L, P, E>(
    nseq: usize,
    selfscore: S,
    orilen: L,
    pair_score: P,
    seqs_equal: E,
    picksize: usize,
) -> Option<PartTreeFromalnResult>
where
    S: Fn(usize) -> i64,
    L: Fn(usize) -> usize,
    P: Fn(usize, usize) -> f64,
    E: Fn(usize, usize) -> bool,
{
    use crate::parttree_pivot::ScoreEntry;
    if nseq <= 1 { return None; }

    // 1) Build initial `entries` with numinseq=0..nseq.
    let mut entries: Vec<ScoreEntry> = (0..nseq).map(|i| {
        ScoreEntry {
            numinseq: i,
            selfscore: selfscore(i),
            orilen: orilen(i),
            score: 0.0,
            points: Vec::new(),
        }
    }).collect();

    // 2) pick_reference: scan for max selfscore (strict `>`, first wins),
    //    swap winner to position 0.
    let mut best = 0usize;
    let mut best_self = entries[0].selfscore;
    for i in 1..nseq {
        if entries[i].selfscore > best_self {
            best_self = entries[i].selfscore;
            best = i;
        }
    }
    if best != 0 { entries.swap(0, best); }

    // 3) Compute scores against entries[0].
    let ref_num = entries[0].numinseq;
    let ref_self = entries[0].selfscore as f64;
    for i in 0..nseq {
        let pair = pair_score(ref_num, entries[i].numinseq);
        let bunbo = ref_self.min(entries[i].selfscore as f64);
        entries[i].score = if bunbo > 0.0 { 1.0 - pair / bunbo } else { 1.0 };
    }
    crate::parttree_pivot::dcompare_sort(&mut entries);

    // 4) Pivot selection (same algorithm as `compute_parttree_order_fromaln`).
    let nin = entries.len();
    let mut picks: Vec<usize> = vec![0];
    let same_seq_at = |a: usize, b: usize| -> bool {
        entries[a].selfscore == entries[b].selfscore
            && entries[a].orilen == entries[b].orilen
            && seqs_equal(entries[a].numinseq, entries[b].numinseq)
    };
    let mut pickkouho: Vec<usize> = (1..nin).collect();
    let mut nkouho = pickkouho.len();
    if nkouho > 0 {
        let picktmp = pickkouho[nkouho - 1];
        nkouho -= 1;
        if !same_seq_at(0, picktmp) { picks.push(picktmp); }
    }
    let mut i_alt = 1;
    while picks.len() < picksize && nkouho > 0 {
        let rn = if i_alt == 1 { i_alt = 0; (nkouho as f64 * 0.5) as usize } else { nkouho - 1 };
        let picktmp = pickkouho[rn];
        nkouho -= 1;
        pickkouho[rn] = pickkouho[nkouho];
        if !picks.iter().any(|&p| same_seq_at(p, picktmp)) {
            picks.push(picktmp);
        }
    }
    picks.sort_unstable();

    // 5) pickmtx via pair_score (same formula as fromaln pipeline).
    let npick = picks.len();
    let mut pickmtx: Vec<Vec<f64>> = (0..npick).map(|i| vec![0.0f64; npick - i]).collect();
    for k in 1..npick {
        pickmtx[0][k] = entries[picks[k]].score;
    }
    for j in 1..npick {
        let pj_self = entries[picks[j]].selfscore as f64;
        let pj_num = entries[picks[j]].numinseq;
        for i in (j + 1)..npick {
            let pair = pair_score(pj_num, entries[picks[i]].numinseq);
            let bunbo = pj_self.min(entries[picks[i]].selfscore as f64);
            let dist = if bunbo > 0.0 { 1.0 - pair / bunbo } else { 1.0 };
            pickmtx[j][i - j] = if dist < 0.0 { 0.0 } else { dist };
        }
    }

    // 6) dfromc[yuko][seq].
    let nyuko = npick;
    let yukos: Vec<usize> = picks.clone();
    let mut dfromc: Vec<Vec<f64>> = vec![vec![0.0f64; nin]; nyuko];
    for j in 0..nin { dfromc[0][j] = entries[j].score; }
    for i in 1..nyuko {
        let yi_num = entries[yukos[i]].numinseq;
        let yi_self = entries[yukos[i]].selfscore as f64;
        for j in 0..nin {
            if j == yukos[i] { dfromc[i][j] = 0.0; continue; }
            let pair = pair_score(yi_num, entries[j].numinseq);
            let bunbo = yi_self.min(entries[j].selfscore as f64);
            let dist = if bunbo > 0.0 { 1.0 - pair / bunbo } else { 1.0 };
            dfromc[i][j] = if dist < 0.0 { 0.0 } else { dist };
        }
    }

    // 7) Yuko assignment.
    let mut outs: Vec<Vec<usize>> = vec![Vec::new(); nyuko];
    for j in 0..nin {
        let mut belongto = 0usize;
        let mut min_d = f64::INFINITY;
        for yi in 0..nyuko {
            if dfromc[yi][j] < min_d {
                min_d = dfromc[yi][j];
                belongto = yi;
            }
        }
        outs[belongto].push(entries[j].numinseq);
    }

    // 8) UPGMA on yukomtx (= pickmtx since every pick survived).
    let mut yuko_dm = DistanceMatrix::new(nyuko);
    for i in 0..nyuko {
        for j in (i + 1)..nyuko {
            yuko_dm.set(i, j, pickmtx[i][j - i]);
        }
    }
    let yuko_topo = musclesupg(&yuko_dm, ClusterMethod::Mix { sueff: 0.1 });

    Some(PartTreeFromalnResult { scores: entries, outs, yuko_topo })
}

/// Given a precomputed parttree pipeline result, build the `--treeout`
/// Newick string (numeric leaves, no branch lengths) — the format C's
/// `splittbfast.c::splitseq_mq` emits via `splittbfast.c:1275-1301`
/// (leaf) and `:2532-2553` (per-merge concat).
pub fn parttree_result_to_newick(result: &PartTreeFromalnResult) -> String {
    let nyuko = result.outs.len();
    let mut parttree: Vec<Option<String>> = (0..nyuko).map(|yi| {
        let members = &result.outs[yi];
        if members.is_empty() {
            None
        } else if members.len() == 1 {
            Some(format!("\n{}\n", members[0] + 1))
        } else {
            let inner: Vec<String> = members.iter().map(|&m| (m + 1).to_string()).collect();
            Some(format!("\n({})\n", inner.join(",")))
        }
    }).collect();

    let mut last_v1 = 0usize;
    for step in &result.yuko_topo.steps {
        let l_norm = c_normalized_subtree(&result.yuko_topo.steps, &step.left);
        let r_norm = c_normalized_subtree(&result.yuko_topo.steps, &step.right);
        let (mut v1, mut v2) = (l_norm[0], r_norm[0]);
        if v1 > v2 { std::mem::swap(&mut v1, &mut v2); }
        let s1 = parttree[v1].take().expect("missing parttree[v1]");
        let s2 = parttree[v2].take().expect("missing parttree[v2]");
        parttree[v1] = Some(format!("({},{})", s1, s2));
        last_v1 = v1;
    }
    let mut out = parttree[last_v1].take().expect("missing root parttree");
    out.push('\n');
    out
}

/// The `--reorder` order of one `splittbfast -Z` (`fromaln=1`) call on
/// `aligned_seqs` (`scripts/mafft:2681`): the leaf emission order of the
/// recursive `splitseq_mq` with `naivepairscore11` distances.
pub fn compute_parttree_order_fromaln(
    aligned_seqs: &[Vec<u8>],
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penalty: f64,
) -> Vec<usize> {
    let nseq = aligned_seqs.len();
    if nseq <= 1 {
        return (0..nseq).collect();
    }
    crate::parttree_splitseq::run_fromaln(aligned_seqs, matrix, amino_map, penalty).order
}

/// PartTree `--treeout` output of one `splittbfast -Z` call on the
/// already-aligned MSA. C overwrites `infile.tree` on each `splittbfast`
/// invocation, so the tree on disk comes from the last (`fromaln=1`) pass.
pub fn compute_parttree_newick_fromaln(
    aligned_seqs: &[Vec<u8>],
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penalty: f64,
) -> String {
    let nseq = aligned_seqs.len();
    if nseq == 0 { return "\n".to_string(); }
    if nseq == 1 { return "\n1\n".to_string(); }
    let mut out = crate::parttree_splitseq::run_fromaln(aligned_seqs, matrix, amino_map, penalty).tree;
    out.push('\n');
    out
}
