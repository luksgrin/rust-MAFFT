/// Pairwise local alignment constraint builder.
///
/// Computes all-vs-all pairwise local alignments and stores the results
/// as a `LocalHomologyTable`, which is then used to guide progressive
/// alignment in L-INS-i and E-INS-i modes.

use rayon::prelude::*;

use mafft_types::{HomologyRegion, LocalHomologyTable};

use crate::dp::GapModel;
use crate::local::local_align;

/// MAFFT's default `fastathreshold` for protein alignments (`scripts/mafft:97`).
/// Multiplies each `region.importance * eff1[gi] * eff2[gj]` contribution
/// before it lands in the per-cell importance matrix.
pub const FASTATHRESHOLD_DEFAULT: f64 = 2.7;

/// Build the per-cell importance matrix `impmtx[i][j]`, mirroring C's
/// `fillimp` (`mltaln9.c:15560-15675`).
///
/// For each pair of group members `(s1, s2)`, walk every homology region
/// stored in `localhom.get(s1, s2)`. The region's `start1`/`end1` are in
/// the original raw-sequence position space; `move_to_seq_pos` translates
/// them to positions in the (gapped) sequences passed in `g1_seqs[gi]` /
/// `g2_seqs[gj]`. Within `[start, end]` we step both sequences in lockstep:
/// at every column where neither has a gap, add
/// `region.importance * eff1[gi] * eff2[gj] * fastathreshold` to
/// `impmtx[k1][k2]`.
///
/// `g1_seqs[gi].len()` must equal `lgth1` (same for group 2). Position
/// translation handles within-group gaps by counting non-gap residues.
pub fn build_imp_matrix(
    localhom: &LocalHomologyTable,
    group1: &[usize],
    group2: &[usize],
    g1_seqs: &[&[u8]],
    g2_seqs: &[&[u8]],
    eff1: &[f64],
    eff2: &[f64],
    lgth1: usize,
    lgth2: usize,
    fastathreshold: f64,
) -> Vec<Vec<f64>> {
    let mut imp = vec![vec![0.0f64; lgth2]; lgth1];
    let effijx = fastathreshold;

    for (gi, &s1) in group1.iter().enumerate() {
        for (gj, &s2) in group2.iter().enumerate() {
            let regions = localhom.get(s1, s2);
            if regions.is_empty() { continue; }
            let effij = eff1[gi] * eff2[gj] * effijx;
            let seq1 = g1_seqs[gi];
            let seq2 = g2_seqs[gj];

            for region in regions {
                let start1 = match move_to_seq_pos(seq1, region.start1 as usize) {
                    Some(p) => p,
                    None => continue,
                };
                let end1 = if region.start1 == region.end1 {
                    start1
                } else {
                    match move_to_seq_pos(seq1, region.end1 as usize) {
                        Some(p) => p,
                        None => continue,
                    }
                };
                let start2 = match move_to_seq_pos(seq2, region.start2 as usize) {
                    Some(p) => p,
                    None => continue,
                };
                let end2 = if region.start2 == region.end2 {
                    start2
                } else {
                    match move_to_seq_pos(seq2, region.end2 as usize) {
                        Some(p) => p,
                        None => continue,
                    }
                };

                let mut k1 = start1;
                let mut k2 = start2;
                while k1 < seq1.len() && k2 < seq2.len() {
                    let c1 = seq1[k1];
                    let c2 = seq2[k2];
                    if c1 != b'-' && c2 != b'-' {
                        if k1 < lgth1 && k2 < lgth2 {
                            imp[k1][k2] += region.importance * effij;
                        }
                        k1 += 1;
                        k2 += 1;
                    } else if c1 != b'-' && c2 == b'-' {
                        k2 += 1;
                    } else if c1 == b'-' && c2 != b'-' {
                        k1 += 1;
                    } else {
                        k1 += 1;
                        k2 += 1;
                    }
                    if k1 > end1 || k2 > end2 { break; }
                }
            }
        }
    }
    imp
}

/// Recompute homology-region `importance` values using C's
/// position-vote algorithm in `calcimportance` (`mltaln9.c:11984`).
///
/// Replaces each region's provisional importance (typically `opt /
/// overlapaa` from `build_local_homology_table`) with a value that
/// rewards regions backed by many "voting" sequences:
///
/// 1. Normalize `ieff[j] = eff[j] / sum(eff[non-empty])`.
/// 2. For each sequence `i`, build a per-position support array
///    `support[pos] = sum_j ieff[j]` over every region in
///    `localhom.get(i, j)` that covers position `pos` (raw residue
///    index in `sequences[i]`).
/// 3. For each region in `localhom.get(i, j)`, set
///    `region.importance = mean(support[start1..=end1]) * region.opt`.
/// 4. Symmetrize: for each pair `(i, j)`, average the importance of
///    the matching regions in `localhom.get(i, j)` and `localhom.get(j, i)`.
///
/// `eff` is the global tree-weight vector (typically from
/// `mafft-tree::sequence_weights`) — match the value C's `tbfast` passes
/// to `calcimportance` after building the tree.
pub fn recompute_importance(
    localhom: &mut LocalHomologyTable,
    sequences: &[&[u8]],
    eff: &[f64],
) {
    let nseq = sequences.len();
    if nseq < 2 || localhom.nseq != nseq { return; }

    let nogaplen: Vec<usize> = sequences.iter()
        .map(|s| s.iter().filter(|&&c| c != b'-').count()).collect();
    let totaleff: f64 = (0..nseq)
        .filter(|&i| nogaplen[i] > 0).map(|i| eff[i]).sum();
    if totaleff <= 0.0 { return; }
    let ieff: Vec<f64> = (0..nseq).map(|i| {
        if nogaplen[i] > 0 { eff[i] / totaleff } else { 0.0 }
    }).collect();
    let nlenmax = nogaplen.iter().copied().max().unwrap_or(0);
    if nlenmax == 0 { return; }

    // Pass 1: per-i position-vote → set importance = mean(support over region) * opt
    let mut support = vec![0.0f64; nlenmax];
    for i in 0..nseq {
        for v in support.iter_mut() { *v = 0.0; }
        for j in 0..nseq {
            if i == j { continue; }
            for region in localhom.get(i, j) {
                let s = (region.start1 as usize).min(nlenmax);
                let e = (region.end1 as usize).min(nlenmax.saturating_sub(1));
                for pos in s..=e {
                    if pos < nlenmax { support[pos] += ieff[j]; }
                }
            }
        }
        for j in 0..nseq {
            if i == j { continue; }
            // Note: C iterates `for tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next`
            // — a linked list. We have a Vec; iterate by index so we can mutate.
            let regions_count = localhom.get(i, j).len();
            for r in 0..regions_count {
                let (s, e) = {
                    let region = &localhom.get(i, j)[r];
                    (
                        (region.start1 as usize).min(nlenmax.saturating_sub(1)),
                        (region.end1 as usize).min(nlenmax.saturating_sub(1)),
                    )
                };
                let mut sum = 0.0f64;
                let mut count = 0usize;
                for pos in s..=e {
                    sum += support[pos];
                    count += 1;
                }
                let mean = if count > 0 { sum / count as f64 } else { 0.0 };
                let region = &mut localhom.get_mut(i, j)[r];
                region.importance = mean * region.opt;
            }
        }
    }

    // Pass 2: symmetrize importance between (i,j) and (j,i)
    for i in 0..nseq.saturating_sub(1) {
        for j in (i + 1)..nseq {
            let n_ij = localhom.get(i, j).len();
            let n_ji = localhom.get(j, i).len();
            let n = n_ij.min(n_ji);
            for r in 0..n {
                let avg = 0.5 * (
                    localhom.get(i, j)[r].importance
                  + localhom.get(j, i)[r].importance
                );
                localhom.get_mut(i, j)[r].importance = avg;
                localhom.get_mut(j, i)[r].importance = avg;
            }
        }
    }
}

/// Translate a raw-sequence residue index (`raw_pos`) into a position in a
/// (possibly gapped) sequence. Mirrors C's `movereg` (`mltaln9.c:15467`):
/// walk forward, count non-gap characters, return the index of the
/// `raw_pos`-th residue. If `raw_pos` exceeds the sequence's residue count
/// (e.g. an exclusive end-marker past the last residue), returns the index
/// just after the last residue, so callers using `<=` end checks still
/// terminate cleanly.
fn move_to_seq_pos(seq: &[u8], raw_pos: usize) -> Option<usize> {
    let target = raw_pos as i64;
    let mut count: i64 = -1;
    let mut last_residue_idx = 0usize;
    let mut found_any = false;
    for (idx, &c) in seq.iter().enumerate() {
        if c != b'-' {
            count += 1;
            last_residue_idx = idx;
            found_any = true;
        }
        if count == target { return Some(idx); }
    }
    if found_any && raw_pos > count as usize {
        Some(last_residue_idx + 1)
    } else {
        None
    }
}

/// Result of one pairwise alignment for parallel collection.
struct PairResult {
    i: usize,
    j: usize,
    distance: f64,
    regions: Vec<HomologyRegion>,
}

/// Which pairwise aligner to drive `build_homology_table` with.
///
/// L-INS-i uses `Local` (Smith–Waterman, mirroring C's `pairlocalalign -L`
/// → `L__align11`). G-INS-i uses `Global` (Needleman–Wunsch, mirroring
/// `pairlocalalign -A` → `G__align11`). E-INS-i uses
/// `GeneralizedAffine` (Smith-Waterman with an extra "skip" gap state,
/// mirroring `pairlocalalign -N` → `genL__align11`). The chaining/`opt`
/// computation downstream is identical — only the alignment differs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairAligner {
    Local,
    Global,
    GeneralizedAffine,
}


/// Build a local homology table from all-vs-all pairwise alignments.
///
/// Defaults to local Smith-Waterman. See `build_homology_table` for the
/// unified entry point that chooses between local and global.
pub fn build_local_homology_table(
    sequences: &[&[u8]],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap: &GapModel,
    score_offset: f64,
) -> (LocalHomologyTable, Vec<Vec<f64>>) {
    build_homology_table(
        sequences, matrix, amino_map, gap, score_offset,
        PairAligner::Local, 0.0,
    )
}

/// Build a homology table using all-vs-all pairwise alignment of the
/// chosen kind. Used by L-INS-i (`PairAligner::Local`), G-INS-i
/// (`PairAligner::Global`), and E-INS-i (`PairAligner::GeneralizedAffine`).
///
/// `op_penalty` is the generalized-affine "skip" open penalty (C's
/// `penalty_OP`), used only when `aligner == GeneralizedAffine`.
pub fn build_homology_table(
    sequences: &[&[u8]],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap: &GapModel,
    score_offset: f64,
    aligner: PairAligner,
    op_penalty: f64,
) -> (LocalHomologyTable, Vec<Vec<f64>>) {
    let nseq = sequences.len();

    // C's `pairlocalalign.c:2590-2596`: selfscore[i] = sum of diagonal
    // substitution-matrix entries for each residue in the sequence,
    // using the same offset-shifted matrix that pairwise alignment uses.
    // Used by `score2dist(pscore, selfscore[i], selfscore[j])`
    // (line 1931) to convert alignment scores to tree-input distances:
    //   bunbo = min(selfscore[i], selfscore[j])
    //   dist = (1 - pscore / bunbo) * 2  (clamped to [0, 2])
    let n_alpha = matrix.len();
    let selfscore: Vec<f64> = sequences.iter().map(|s| {
        let mut sum = 0.0f64;
        for &c in *s {
            let i = amino_map[c as usize] as usize;
            if i < n_alpha {
                sum += matrix[i][i] as f64;
            }
        }
        sum
    }).collect();

    // Generate all (i, j) pairs with i < j
    let pairs: Vec<(usize, usize)> = (0..nseq)
        .flat_map(|i| ((i + 1)..nseq).map(move |j| (i, j)))
        .collect();

    // Compute all pairwise alignments in parallel
    let results: Vec<PairResult> = pairs
        .par_iter()
        .map(|&(i, j)| {
            // Branch at the alignment call to keep the rest of the
            // chaining logic shared. For Local, we have a true offset
            // into both sequences. For Global, both offsets are 0.
            let (alignment, offset1, offset2) = match aligner {
                PairAligner::Local => {
                    let r = local_align(
                        sequences[i], sequences[j],
                        matrix, amino_map, gap, score_offset,
                    );
                    (r.alignment, r.offset1, r.offset2)
                }
                PairAligner::Global => {
                    // Match C's `pairlocalalign -A` (`pairlocalalign.c:2196`)
                    // which calls `G__align11` with `outgap` controlling
                    // terminal-gap treatment. The script for G-INS-i
                    // doesn't pass `-O`, so `outgap` defaults to 1 →
                    // both head and tail gaps are penalized.
                    let r = crate::global::global_align(
                        sequences[i], sequences[j],
                        matrix, amino_map, gap,
                        true, true,
                    );
                    (r, 0, 0)
                }
                PairAligner::GeneralizedAffine => {
                    // C's `pairlocalalign -N` (`pairlocalalign.c:2233`) →
                    // `genL__align11`: max-so-far Smith-Waterman with an
                    // extra "skip" gap state (penalty_OP, no extension).
                    use crate::genaffine::{genaffine_local_align, GenAffineGapModel};
                    let gen_gap = GenAffineGapModel {
                        affine: gap.clone(),
                        open_generalized: op_penalty,
                    };
                    let r = genaffine_local_align(
                        sequences[i], sequences[j],
                        matrix, amino_map, &gen_gap, score_offset,
                    );
                    (r.alignment, r.offset1, r.offset2)
                }
            };

            let score = alignment.score;
            if score <= 0.0 {
                return PairResult { i, j, distance: 2.0, regions: Vec::new() };
            }

            // C's `score2dist` (`pairlocalalign.c:1931-1944`):
            //   bunbo = min(selfscore[i], selfscore[j])
            //   dist = bunbo == 0           ? 2.0
            //        : bunbo < pscore       ? 0.0
            //        : (1 - pscore / bunbo) * 2
            let bunbo = selfscore[i].min(selfscore[j]);
            let d = if bunbo == 0.0 {
                2.0
            } else if bunbo < score {
                0.0
            } else {
                (1.0 - score / bunbo) * 2.0
            };

            // Port of C's `putlocalhom2` (`io.c:723`): split the alignment
            // into maximal gap-free regions. Whenever a gap appears in
            // either aligned sequence we close the current region; a
            // residue-residue column starts a new one. Each region
            // records its (start1, end1, start2, end2) span and its
            // contribution `iscore` to `sumoverlap` / `isumscore`.
            //
            // Since L-INS-i runs with `divpairscore = 0` (no `-y`), every
            // region in the chain gets the same combined opt. C's
            // `putlocalhom2` stores `isumscore * 5.8 / (600 * sumoverlap)`
            // (normalized for hat3 file output), then `tbfast.c:2202`
            // immediately rescales it back via `opt * 600 / 5.8` before
            // calling `calcimportance_half`. We bypass that round-trip and
            // store the post-scale value `isumscore / sumoverlap`
            // directly, which matches the value C's calcimportance/fillimp
            // actually consume.
            let n_alpha = matrix.len();
            let a1 = &alignment.seq1;
            let a2 = &alignment.seq2;
            let mut regions: Vec<HomologyRegion> = Vec::new();
            let mut isumscore: f64 = 0.0;
            let mut sumoverlap: i32 = 0;
            let mut pos1 = offset1 as i32;
            let mut pos2 = offset2 as i32;
            let mut st = false;
            let mut start1 = 0i32;
            let mut start2 = 0i32;
            let mut iscore: f64 = 0.0;
            for k in 0..a1.len() {
                let c1 = a1[k];
                let c2 = a2[k];
                let g1 = c1 == b'-';
                let g2 = c2 == b'-';
                if st && (g1 || g2) {
                    let end1 = pos1 - 1;
                    let end2 = pos2 - 1;
                    regions.push(HomologyRegion {
                        start1, end1, start2, end2,
                        opt: 0.0,                       // filled below
                        overlapaa: end2 - start2 + 1,
                        korh: b'h',
                        ..Default::default()
                    });
                    isumscore += iscore;
                    sumoverlap += end2 - start2 + 1;
                    iscore = 0.0;
                    st = false;
                } else if !g1 && !g2 {
                    if !st {
                        start1 = pos1;
                        start2 = pos2;
                        st = true;
                    }
                    let i1 = amino_map[c1 as usize] as usize;
                    let i2 = amino_map[c2 as usize] as usize;
                    if i1 < n_alpha && i2 < n_alpha {
                        iscore += matrix[i1][i2] as f64;
                    }
                }
                if !g1 { pos1 += 1; }
                if !g2 { pos2 += 1; }
            }
            // Close trailing region if alignment ends inside a match span.
            if st {
                let end1 = pos1 - 1;
                let end2 = pos2 - 1;
                regions.push(HomologyRegion {
                    start1, end1, start2, end2,
                    opt: 0.0,
                    overlapaa: end2 - start2 + 1,
                    korh: b'h',
                    ..Default::default()
                });
                isumscore += iscore;
                sumoverlap += end2 - start2 + 1;
            }

            // !divpairscore branch (`io.c:855-866` + `tbfast.c:2202` rescale):
            // all regions share a single combined opt and overlapaa.
            let opt = if sumoverlap > 0 {
                isumscore / sumoverlap as f64
            } else { 0.0 };
            let provisional_importance =
                if sumoverlap > 0 { opt / sumoverlap as f64 } else { 0.0 };
            for r in regions.iter_mut() {
                r.opt = opt;
                r.overlapaa = sumoverlap;
                r.importance = provisional_importance;
            }

            PairResult { i, j, distance: d, regions }
        })
        .collect();

    // Apply results to table and distance matrix (sequential)
    let mut table = LocalHomologyTable::new(nseq);
    let mut dist = vec![vec![0.0f64; nseq]; nseq];

    for r in results {
        dist[r.i][r.j] = r.distance;
        dist[r.j][r.i] = r.distance;

        for region in &r.regions {
            table.push(r.i, r.j, region.clone());
            table.push(r.j, r.i, HomologyRegion {
                start1: region.start2,
                end1: region.end2,
                start2: region.start1,
                end2: region.end1,
                ..region.clone()
            });
        }
    }

    (table, dist)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_setup() -> (Vec<Vec<i32>>, [u8; 256]) {
        let mut mtx = vec![vec![-100i32; 5]; 5];
        for i in 0..4 { mtx[i][i] = 100; }
        let mut map = [0xFFu8; 256];
        map[b'A' as usize] = 0; map[b'C' as usize] = 1;
        map[b'G' as usize] = 2; map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;
        (mtx, map)
    }

    #[test]
    fn builds_table_for_similar_sequences() {
        let (mtx, map) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGTACGT"];
        let gap = GapModel::new(-200.0, -10.0);

        let (table, dist) = build_local_homology_table(&seqs, &mtx, &map, &gap, 0.0);

        let regions_01 = table.get(0, 1);
        assert!(!regions_01.is_empty(), "should find homology between identical seqs");
        assert!(dist[0][1] < 1e-6, "identical seqs should have distance ~0, got {}", dist[0][1]);
    }

    #[test]
    fn reciprocal_entries() {
        let (mtx, map) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGTACGT"];
        let gap = GapModel::new(-200.0, -10.0);

        let (table, _) = build_local_homology_table(&seqs, &mtx, &map, &gap, 0.0);

        let fwd = table.get(0, 1);
        let rev = table.get(1, 0);
        assert_eq!(fwd.len(), rev.len());

        if !fwd.is_empty() {
            assert_eq!(fwd[0].start1, rev[0].start2);
            assert_eq!(fwd[0].start2, rev[0].start1);
        }
    }
}
