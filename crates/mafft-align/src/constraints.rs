/// Pairwise local alignment constraint builder.
///
/// Computes all-vs-all pairwise local alignments and stores the results
/// as a `LocalHomologyTable`, which is then used to guide progressive
/// alignment in L-INS-i and E-INS-i modes.

use rayon::prelude::*;

use mafft_types::{HomologyRegion, LocalHomologyTable};

use crate::dp::GapModel;
use crate::local::local_align;

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
