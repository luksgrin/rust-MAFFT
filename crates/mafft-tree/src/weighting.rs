/// Sequence weighting from tree topology.
///
/// Ports the C `counteff_simple_double()` from mltaln9.c.
/// Down-weights redundant sequences based on their position in the guide tree.

use crate::topology::Topology;

/// Small constant added to all weights to prevent zero weights.
const GETA: f64 = 0.001;

/// Compute sequence weights from a guide tree topology.
///
/// Sequences that are more isolated in the tree (longer branches) get
/// higher weights, while closely related sequences get lower weights.
/// This improves alignment quality for datasets with uneven sampling.
///
/// Returns a weight vector of length `nseq`, normalized to sum to 1.0.
pub fn sequence_weights(topo: &Topology) -> Vec<f64> {
    let n = topo.nseq;
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![1.0];
    }

    let mut rootnode = vec![0.0f64; n];
    let mut eff = vec![1.0f64; n];

    for step in &topo.steps {
        // Left cluster: accumulate branch length weighted by current efficiency
        for &s in &step.left {
            rootnode[s] += step.left_length * eff[s];
            eff[s] *= 0.5;
        }

        // Right cluster
        for &s in &step.right {
            rootnode[s] += step.right_length * eff[s];
            eff[s] *= 0.5;
        }
    }

    // Add small constant and normalize
    for w in &mut rootnode {
        *w += GETA;
    }

    let total: f64 = rootnode.iter().sum();
    if total > 0.0 {
        for w in &mut rootnode {
            *w /= total;
        }
    }

    rootnode
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::JoinStep;

    #[test]
    fn weights_sum_to_one() {
        let mut topo = Topology::new(3);
        topo.steps.push(JoinStep {
            left: vec![0],
            right: vec![1],
            left_length: 0.1,
            right_length: 0.2,
        });
        topo.steps.push(JoinStep {
            left: vec![0, 1],
            right: vec![2],
            left_length: 0.3,
            right_length: 0.5,
        });

        let w = sequence_weights(&topo);
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "weights sum to {sum}");
    }

    #[test]
    fn isolated_sequence_gets_higher_weight() {
        // Seqs 0,1 are close; seq 2 is distant
        let mut topo = Topology::new(3);
        topo.steps.push(JoinStep {
            left: vec![0],
            right: vec![1],
            left_length: 0.05,
            right_length: 0.05,
        });
        topo.steps.push(JoinStep {
            left: vec![0, 1],
            right: vec![2],
            left_length: 0.1,
            right_length: 0.9,
        });

        let w = sequence_weights(&topo);
        // Seq 2 (isolated) should have higher weight than 0 or 1
        assert!(w[2] > w[0], "isolated seq should have higher weight");
        assert!(w[2] > w[1], "isolated seq should have higher weight");
    }

    #[test]
    fn single_sequence() {
        let w = sequence_weights(&Topology::new(1));
        assert_eq!(w, vec![1.0]);
    }
}
