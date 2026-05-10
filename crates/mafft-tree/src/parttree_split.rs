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
