/// Pairwise local alignment constraint builder.
///
/// Ports the C `pairlocalalign.c`.
///
/// Computes all-vs-all pairwise local alignments and stores the results
/// as a `LocalHomologyTable`, which is then used to guide progressive
/// alignment in L-INS-i and E-INS-i modes.

use mafft_types::{HomologyRegion, LocalHomologyTable};

use crate::dp::GapModel;
use crate::local::local_align;

/// Build a local homology table from all-vs-all pairwise local alignments.
///
/// For each pair (i, j) with i < j, performs a local alignment and stores
/// the resulting high-scoring regions. Reciprocal entries (j, i) are created
/// with swapped coordinates.
///
/// Returns the table and a pairwise distance matrix.
pub fn build_local_homology_table(
    sequences: &[&[u8]],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap: &GapModel,
    score_offset: f64,
) -> (LocalHomologyTable, Vec<Vec<f64>>) {
    let nseq = sequences.len();
    let mut table = LocalHomologyTable::new(nseq);
    let mut dist = vec![vec![0.0f64; nseq]; nseq];

    for i in 0..nseq {
        for j in (i + 1)..nseq {
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
                dist[i][j] = 2.0;
                dist[j][i] = 2.0;
                continue;
            }

            // Compute distance: use identity from the alignment
            let identity = result.alignment.identity();
            let d = (1.0 - identity).clamp(0.0, 2.0);
            dist[i][j] = d;
            dist[j][i] = d;

            // Store local homology region
            let aln_len = result.alignment.len();
            if aln_len > 0 {
                let region = HomologyRegion {
                    start1: result.offset1 as i32,
                    end1: (result.offset1 + aln_len) as i32,
                    start2: result.offset2 as i32,
                    end2: (result.offset2 + aln_len) as i32,
                    opt: score,
                    overlapaa: aln_len as i32,
                    importance: score,
                    korh: b'h',
                    ..Default::default()
                };
                table.push(i, j, region.clone());

                // Reciprocal entry
                table.push(j, i, HomologyRegion {
                    start1: region.start2,
                    end1: region.end2,
                    start2: region.start1,
                    end2: region.end1,
                    ..region
                });
            }
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
        map[b'A' as usize] = 0;
        map[b'C' as usize] = 1;
        map[b'G' as usize] = 2;
        map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;
        (mtx, map)
    }

    #[test]
    fn builds_table_for_similar_sequences() {
        let (mtx, map) = simple_setup();
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGTACGT"];
        let gap = GapModel::new(-200.0, -10.0);

        let (table, dist) = build_local_homology_table(&seqs, &mtx, &map, &gap, 0.0);

        // Identical sequences should have homology regions
        let regions_01 = table.get(0, 1);
        assert!(!regions_01.is_empty(), "should find homology between identical seqs");

        // Distance between identical seqs should be 0
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
            // Reciprocal should have swapped coordinates
            assert_eq!(fwd[0].start1, rev[0].start2);
            assert_eq!(fwd[0].start2, rev[0].start1);
        }
    }
}
