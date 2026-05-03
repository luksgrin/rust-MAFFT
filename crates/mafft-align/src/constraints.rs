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
    region: Option<HomologyRegion>,
}

/// Build a local homology table from all-vs-all pairwise local alignments.
///
/// The pairwise alignments are computed in parallel using rayon.
/// Results are collected and applied to the table sequentially.
pub fn build_local_homology_table(
    sequences: &[&[u8]],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap: &GapModel,
    score_offset: f64,
) -> (LocalHomologyTable, Vec<Vec<f64>>) {
    let nseq = sequences.len();

    // Generate all (i, j) pairs with i < j
    let pairs: Vec<(usize, usize)> = (0..nseq)
        .flat_map(|i| ((i + 1)..nseq).map(move |j| (i, j)))
        .collect();

    // Compute all pairwise alignments in parallel
    let results: Vec<PairResult> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let result = local_align(
                sequences[i],
                sequences[j],
                matrix,
                amino_map,
                gap,
                score_offset,
            );

            let score = result.alignment.score;
            if score <= 0.0 {
                return PairResult { i, j, distance: 2.0, region: None };
            }

            let identity = result.alignment.identity();
            let d = (1.0 - identity).clamp(0.0, 2.0);

            let aln_len = result.alignment.len();

            // C's `pairlocalalign` recomputes `opt` from the alignment
            // by walking the matched residues and summing
            // `n_dis[c1][c2]` cells, then rescales:
            //   opt = iscore * 5.8 / (600 * sumoverlap)
            // (`pairlocalalign.c:201` for divpairscore=1, line 222 for
            // divpairscore=0.) The `score` we get from `local_align` is
            // the full DP score (matches plus gap penalties), so we
            // recompute `iscore` from the alignment to match C exactly.
            let mut iscore: f64 = 0.0;
            let n_alpha = matrix.len();
            for k in 0..result.alignment.seq1.len() {
                let c1 = result.alignment.seq1[k];
                let c2 = result.alignment.seq2[k];
                if c1 != b'-' && c2 != b'-' {
                    let i1 = amino_map[c1 as usize] as usize;
                    let i2 = amino_map[c2 as usize] as usize;
                    if i1 < n_alpha && i2 < n_alpha {
                        iscore += matrix[i1][i2] as f64;
                    }
                }
            }
            let opt = if aln_len > 0 {
                iscore * 5.8 / (600.0 * aln_len as f64)
            } else { 0.0 };

            // Provisional importance = opt / overlapaa (`dontcalcimportance`,
            // mltaln9.c:11472). `calcimportance_half` (mltaln9.c:11756)
            // overwrites this with `mean(support) * opt` after the tree is
            // built — done in `recompute_importance` if enabled.
            let importance = if aln_len > 0 { opt / aln_len as f64 } else { 0.0 };
            let region = if aln_len > 0 {
                Some(HomologyRegion {
                    start1: result.offset1 as i32,
                    end1: (result.offset1 + aln_len) as i32,
                    start2: result.offset2 as i32,
                    end2: (result.offset2 + aln_len) as i32,
                    opt,
                    overlapaa: aln_len as i32,
                    importance,
                    korh: b'h',
                    ..Default::default()
                })
            } else {
                None
            };

            PairResult { i, j, distance: d, region }
        })
        .collect();

    // Apply results to table and distance matrix (sequential)
    let mut table = LocalHomologyTable::new(nseq);
    let mut dist = vec![vec![0.0f64; nseq]; nseq];

    for r in results {
        dist[r.i][r.j] = r.distance;
        dist[r.j][r.i] = r.distance;

        if let Some(region) = r.region {
            table.push(r.i, r.j, region.clone());
            table.push(r.j, r.i, HomologyRegion {
                start1: region.start2,
                end1: region.end2,
                start2: region.start1,
                end2: region.end1,
                ..region
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
