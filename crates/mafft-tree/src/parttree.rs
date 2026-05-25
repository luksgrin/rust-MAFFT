/// PartTree: divide-and-conquer guide tree for large datasets (10K+ sequences).
///
/// Ports the C `splitseq_mq()` from `splittbfast.c`.
///
/// Algorithm:
/// 1. Select diverse pivot sequences
/// 2. Compute distances from pivots to all sequences
/// 3. Assign each sequence to its nearest pivot
/// 4. Recursively partition until groups are small enough
/// 5. Build UPGMA trees within small groups
/// 6. Combine partitions into a full topology

use crate::distance::{DistanceMatrix, ktuple_distance};
use crate::musclesupg::{musclesupg, ClusterMethod};
use crate::topology::{JoinStep, Topology};

/// Parameters for PartTree construction.
#[derive(Debug, Clone)]
pub struct PartTreeParams {
    /// Maximum group size before switching to full UPGMA (default: 150).
    pub group_size: usize,
    /// Number of pivot sequences per partition (default: 50).
    pub pick_size: usize,
    /// Use DP alignment for distances instead of k-tuple (--dpparttree).
    pub use_dp: bool,
}

impl Default for PartTreeParams {
    fn default() -> Self {
        Self {
            group_size: 150,
            pick_size: 50,
            use_dp: false,
        }
    }
}

/// Build a guide tree using the PartTree algorithm.
///
/// For datasets with fewer than `group_size` sequences, falls back to
/// standard UPGMA. For larger datasets, recursively partitions into
/// groups and builds trees bottom-up.
pub fn parttree(
    sequences: &[Vec<u8>],
    params: &PartTreeParams,
) -> Topology {
    let n = sequences.len();
    if n <= 1 {
        return Topology::new(n);
    }

    let indices: Vec<usize> = (0..n).collect();

    if n <= params.group_size {
        // Small enough for full UPGMA
        return build_upgma_for_group(sequences, &indices);
    }

    // Recursive partitioning
    let mut topo = Topology::new(n);
    partition_recursive(sequences, &indices, params, &mut topo);
    topo
}

/// Recursively partition sequences and build tree.
fn partition_recursive(
    sequences: &[Vec<u8>],
    indices: &[usize],
    params: &PartTreeParams,
    topo: &mut Topology,
) {
    let n = indices.len();

    if n <= 1 {
        return;
    }

    if n <= params.group_size {
        // Base case: build UPGMA tree for this group
        let sub_topo = build_upgma_for_group(sequences, indices);
        // Append sub-topology steps to main topology
        for step in &sub_topo.steps {
            topo.steps.push(step.clone());
        }
        return;
    }

    // Select pivots
    let npick = params.pick_size.min(n);
    let pivots = select_pivots(sequences, indices, npick);

    if pivots.len() < 2 {
        // Can't partition further, just build UPGMA
        let sub_topo = build_upgma_for_group(sequences, indices);
        for step in &sub_topo.steps {
            topo.steps.push(step.clone());
        }
        return;
    }

    // Compute distances from first two pivots to all sequences
    let pivot0 = pivots[0];
    let pivot1 = pivots[1];

    // Assign each sequence to nearest pivot (binary partition)
    let mut group0 = Vec::new();
    let mut group1 = Vec::new();

    for &idx in indices {
        let d0 = ktuple_distance(&sequences[idx], &sequences[pivot0], 6);
        let d1 = ktuple_distance(&sequences[idx], &sequences[pivot1], 6);
        if d0 <= d1 {
            group0.push(idx);
        } else {
            group1.push(idx);
        }
    }

    // Handle degenerate cases
    if group0.is_empty() || group1.is_empty() {
        let sub_topo = build_upgma_for_group(sequences, indices);
        for step in &sub_topo.steps {
            topo.steps.push(step.clone());
        }
        return;
    }

    // Recurse on both groups
    partition_recursive(sequences, &group0, params, topo);
    partition_recursive(sequences, &group1, params, topo);

    // Merge the two groups at this level
    topo.steps.push(JoinStep {
        left: group0,
        right: group1,
        left_length: 0.5,
        right_length: 0.5,
    });
}

/// Select diverse pivot sequences.
///
/// First pivot: longest sequence. Second pivot: most distant from first.
/// Remaining pivots: maximize minimum distance to already-selected pivots.
fn select_pivots(
    sequences: &[Vec<u8>],
    indices: &[usize],
    npick: usize,
) -> Vec<usize> {
    let n = indices.len();
    let npick = npick.min(n);

    if npick == 0 {
        return Vec::new();
    }

    // First pivot: longest sequence
    let mut pivots = Vec::with_capacity(npick);
    let first = *indices.iter()
        .max_by_key(|&&i| sequences[i].len())
        .unwrap();
    pivots.push(first);

    if npick == 1 {
        return pivots;
    }

    // Compute distances from first pivot to all
    let mut dist_to_nearest_pivot: Vec<f64> = indices.iter()
        .map(|&i| ktuple_distance(&sequences[i], &sequences[first], 6))
        .collect();

    // Greedily select pivots maximizing minimum distance to existing pivots
    for _ in 1..npick {
        // Find sequence with maximum distance to nearest pivot
        let (best_local_idx, _) = dist_to_nearest_pivot.iter()
            .enumerate()
            .filter(|(li, _)| !pivots.contains(&indices[*li]))
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, &0.0));

        let best = indices[best_local_idx];
        pivots.push(best);

        // Update distances: for each sequence, take min distance to any pivot
        for (li, &idx) in indices.iter().enumerate() {
            let d = ktuple_distance(&sequences[idx], &sequences[best], 6);
            if d < dist_to_nearest_pivot[li] {
                dist_to_nearest_pivot[li] = d;
            }
        }
    }

    pivots
}

/// Build a full UPGMA tree for a subset of sequences.
fn build_upgma_for_group(
    sequences: &[Vec<u8>],
    indices: &[usize],
) -> Topology {
    let n = indices.len();
    if n <= 1 {
        return Topology::new(n);
    }

    let mut dm = DistanceMatrix::new(n);
    for i in 0..n {
        for j in (i + 1)..n {
            let d = ktuple_distance(&sequences[indices[i]], &sequences[indices[j]], 6);
            dm.set(i, j, d);
        }
    }

    let local_topo = musclesupg(&dm, ClusterMethod::default());

    // Remap local indices (0..n) to global indices
    let mut remapped = Topology::new(n);
    for step in &local_topo.steps {
        remapped.steps.push(JoinStep {
            left: step.left.iter().map(|&li| indices[li]).collect(),
            right: step.right.iter().map(|&li| indices[li]).collect(),
            left_length: step.left_length,
            right_length: step.right_length,
        });
    }
    remapped
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_seqs(n: usize) -> Vec<Vec<u8>> {
        // Generate n protein-like sequences with varying similarity
        let base = b"ACDEFGHIKLMNPQRSTVWY";
        (0..n).map(|i| {
            let mut seq = Vec::with_capacity(100);
            for j in 0..100 {
                seq.push(base[(j + i * 3) % base.len()]);
            }
            seq
        }).collect()
    }

    #[test]
    fn parttree_small_falls_back_to_upgma() {
        let seqs = make_seqs(10);
        let params = PartTreeParams { group_size: 150, ..Default::default() };
        let topo = parttree(&seqs, &params);
        assert_eq!(topo.nseq, 10);
        // Should have 9 merge steps
        assert_eq!(topo.steps.len(), 9);
        // All sequences should be in the last step
        let last = topo.steps.last().unwrap();
        let mut all: Vec<usize> = last.left.iter().chain(last.right.iter()).copied().collect();
        all.sort();
        assert_eq!(all, (0..10).collect::<Vec<_>>());
    }

    #[test]
    fn parttree_large_partitions() {
        let seqs = make_seqs(200);
        let params = PartTreeParams { group_size: 50, pick_size: 10, use_dp: false };
        let topo = parttree(&seqs, &params);
        assert_eq!(topo.nseq, 200);
        // Should have 199 merge steps
        assert_eq!(topo.steps.len(), 199);
        // All sequences should appear
        let last = topo.steps.last().unwrap();
        let mut all: Vec<usize> = last.left.iter().chain(last.right.iter()).copied().collect();
        all.sort();
        assert_eq!(all, (0..200).collect::<Vec<_>>());
    }

    #[test]
    fn parttree_single_sequence() {
        let seqs = make_seqs(1);
        let topo = parttree(&seqs, &PartTreeParams::default());
        assert_eq!(topo.nseq, 1);
        assert!(topo.steps.is_empty());
    }

    #[test]
    fn parttree_two_sequences() {
        let seqs = make_seqs(2);
        let topo = parttree(&seqs, &PartTreeParams::default());
        assert_eq!(topo.nseq, 2);
        assert_eq!(topo.steps.len(), 1);
    }
}
