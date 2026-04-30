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
            let region = if aln_len > 0 {
                Some(HomologyRegion {
                    start1: result.offset1 as i32,
                    end1: (result.offset1 + aln_len) as i32,
                    start2: result.offset2 as i32,
                    end2: (result.offset2 + aln_len) as i32,
                    opt: score,
                    overlapaa: aln_len as i32,
                    importance: score,
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
