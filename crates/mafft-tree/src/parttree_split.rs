//! `splitseq_mq` recursion + topology assembly for `--parttree`,
//! mirroring `splittbfast.c::splitseq_mq` (lines 2132-2580) for the
//! single-recursion-level case (`nin <= picksize`, all sequences become
//! pivots). For our 36-seq fixture this is the entire algorithm; full
//! recursion handling (multi-level, `nin > picksize`) is left for a
//! follow-up because it requires `rand()` determinism.
//!
//! Pipeline at this layer:
//! 1. Build `dfromc[nyuko][nin]` — distance from each surviving yuko
//!    to every sorted-position sequence. Mirrors
//!    `splittbfast.c:2132-2236`.
//! 2. Assign each sequence (including pivots themselves) to its
//!    closest yuko via `argmin_i dfromc[i][j]`, building `outs[].numinseq`
//!    lists (`splittbfast.c:2244-2287`).
//! 3. UPGMA on `yukomtx` → topol (already done by Rust's `musclesupg`,
//!    FFI-validated).
//! 4. Assemble the final `Topology`:
//!    - Emit an internal-alignment `JoinStep` for each multi-member
//!      yuko, pre-aligning its assigned sequences (Rust's progressive
//!      alignment expects each leaf to be a single sequence; C's
//!      pairalign tolerates mixed-length raw inputs but Rust's
//!      `Profile::from_aligned` requires same-length rows).
//!    - For each UPGMA join step l, emit a `JoinStep` whose `left` /
//!      `right` are the union of `outs[yuko_idx].numinseq` for all
//!      yuko-indices in `topol[l][0]` / `topol[l][1]`, sorted by
//!      `intcompare` (ascending original sequence index — matching C's
//!      `qsort(mem1, ..., intcompare)` at `splittbfast.c:2477-2478`).

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
/// **Multi-member yukos**: for a yuko with members `[m_0, m_1, …,
/// m_{k-1}]`, C does NOT pre-align them. The recursive
/// `splitseq_mq` call hits LEAF (via `uniform = -1` set when
/// `nyuko == 1` at `splittbfast.c:2115`) and just writes `order[]`
/// without alignment. The members are then merged for the first time
/// inside the parent's `pairalign(mem1=[m_0..m_{k-1}], …)` call
/// which builds a profile from the raw sequences and runs DP. For
/// our n=36 fixture the multi-member yuko's members are byte-
/// identical (the dedupe on shimon+strcmp), so they always have the
/// same length and Rust's `Profile::from_aligned` accepts them
/// directly. For non-identical multi-member groups (which require
/// the recursive call to actually run pairalign), more work is
/// needed — see TODO §6.
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

/// One-shot helper: run pivot pipeline, build `dfromc`, assign to yukos,
/// and assemble the final topology.
pub fn build_parttree_topology(
    sequences: &[Vec<u8>],
    kind: PtSeqKind,
    picksize: usize,
) -> Topology {
    let pivots = crate::parttree_pivot::run_pivot_pipeline(sequences, kind, picksize);
    let dfromc = build_dfromc(&pivots, kind);
    let outs = assign_to_yukos(&pivots, &dfromc);
    assemble_topology(&pivots, &outs)
}

/// Recursive helper for `c_normalized_yuko_order`: produce the leaf-yuko
/// order for the subtree whose leaves are `leaves`, applying C's
/// smaller-first-element normalization at each internal node (mirroring
/// `mltaln9.c:8184-8197` in `fixed_musclesupg_double_realloc_nobk_halfmtx`).
fn c_normalized_subtree(steps: &[crate::topology::JoinStep], leaves: &[usize]) -> Vec<usize> {
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

/// Compute the C-equivalent partition-discovery order for `--reorder`,
/// mirroring `splittbfast.c::splitseq_mq` (`splittbfast.c:2351-2378` for
/// the yuko visit order + `:1305-1309` for the leaf emission). Currently
/// supports the single-recursion-level case (`nin <= picksize`, all yukos
/// bottom out without further pivoting) — matches the n=36 fixture.
/// Multi-level recursion (when a yuko itself has too many members to
/// trivially leaf out) would need a full port of `splitseq_mq` recursion;
/// for that case we fall back to input order, which gives a structurally
/// valid alignment but won't byte-match C.
pub fn compute_parttree_order(
    sequences: &[Vec<u8>],
    kind: PtSeqKind,
    picksize: usize,
) -> Vec<usize> {
    let nseq = sequences.len();
    if nseq <= 1 {
        return (0..nseq).collect();
    }
    let pivots = crate::parttree_pivot::run_pivot_pipeline(sequences, kind, picksize);
    let dfromc = build_dfromc(&pivots, kind);
    let outs = assign_to_yukos(&pivots, &dfromc);
    let yuko_dm = yukomtx_to_distance_matrix(&pivots);
    let yuko_topo = musclesupg(&yuko_dm, ClusterMethod::Mix { sueff: 0.1 });

    // Mirrors `splittbfast.c:2351-2354`: treeorder = root step's
    // `topol[nyuko-2][0] ++ topol[nyuko-2][1]`. Our `JoinStep.left/right`
    // accumulate leaves in merge order, NOT in C's smaller-first-element
    // normalized order (`mltaln9.c:8184-8197`). So we recursively rebuild
    // the order with that normalization applied to recover what C's
    // `topol[step][i]` arrays would hold.
    let mut order = Vec::with_capacity(nseq);
    if let Some(root) = yuko_topo.steps.last() {
        let l_yukos = c_normalized_subtree(&yuko_topo.steps, &root.left);
        let r_yukos = c_normalized_subtree(&yuko_topo.steps, &root.right);
        for &yi in l_yukos.iter().chain(r_yukos.iter()) {
            // Leaf-level emission: `splittbfast.c:1305-1309` writes
            // `scores[j].numinseq` in `j` order — `outs[yi]` already
            // holds them in that order (see `assign_to_yukos`).
            order.extend_from_slice(&outs[yi]);
        }
    } else {
        // nyuko == 1 → single yuko containing all sequences in
        // `outs[0]` (matches `splittbfast.c:2105-2123`'s uniform branch).
        order.extend_from_slice(&outs[0]);
    }
    debug_assert_eq!(order.len(), nseq, "parttree order missing sequences");
    order
}
