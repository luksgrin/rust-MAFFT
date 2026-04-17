/// Sequence weighting from tree topology.
///
/// Ports the C `counteff_simple_double()` from mltaln9.c (global weights)
/// and `weightFromABranch()` from treeOperation.c (per-branch weights).

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

// ---------------------------------------------------------------------------
// Per-branch weighting (C's weightFromABranch / calcBranchWeight)
// ---------------------------------------------------------------------------

/// Unrooted tree node for per-branch weight computation.
/// Each internal node has 3 neighbors; leaves have 1 neighbor.
struct WNode {
    /// Indices of neighbor nodes (up to 3). -1 = no neighbor.
    children: [i32; 3],
    /// Branch lengths to each neighbor.
    length: [f64; 3],
    /// Pre-computed branch weights (from calcW) for each neighbor direction.
    branch_weight: [f64; 3],
    /// Sequence members reachable in each direction.
    members: [Vec<i32>; 3],
}

/// Pre-computed per-branch weight structure.
/// Build once from a Topology, then call `weights_for_branch` per branch.
pub struct BranchWeights {
    nodes: Vec<WNode>,
    nseq: usize,
}

impl BranchWeights {
    /// Build the per-branch weight structure from a topology.
    /// Ports C's `treeCnv` + `calcBranchWeight`.
    pub fn new(topo: &Topology) -> Self {
        let nseq = topo.nseq;
        if nseq <= 2 {
            return Self { nodes: Vec::new(), nseq };
        }

        // C allocates 2*nseq nodes (calloc(locnjob*2, sizeof(Node))).
        // Indices: 0..nseq-2 for internal nodes (merge steps),
        // nseq..2*nseq-1 for leaf nodes.
        let total = 2 * nseq;
        let mut nodes: Vec<WNode> = (0..total).map(|_| WNode {
            children: [-1; 3],
            length: [0.0; 3],
            branch_weight: [1.0; 3],
            members: [Vec::new(), Vec::new(), Vec::new()],
        }).collect();

        // Assign leaf members: leaf node i (at index nseq+i) represents sequence i.
        for i in 0..nseq {
            nodes[nseq + i].members[0] = vec![i as i32, -1];
        }

        // Build tree from topology steps.
        // Each step k merges left and right clusters, creating internal node k.
        // We link each leaf/subtree to its parent internal node.
        // Track which internal node each sequence's cluster belongs to.
        let mut parent_of_leaf: Vec<(usize, usize)> = vec![(usize::MAX, 0); nseq]; // (step_idx, 0=left/1=right)

        for (k, step) in topo.steps.iter().enumerate() {
            // Internal node k has children from step.left and step.right.
            // members[0] = left cluster members, members[1] = right cluster members.
            nodes[k].members[0] = step.left.iter().map(|&s| s as i32).chain(std::iter::once(-1)).collect();
            nodes[k].members[1] = step.right.iter().map(|&s| s as i32).chain(std::iter::once(-1)).collect();

            // Find the child nodes for left side.
            // If left is a single sequence, child is the leaf node.
            // If left is a merged cluster from a previous step, child is that step's node.
            let left_child = if step.left.len() == 1 {
                nseq + step.left[0] // leaf node
            } else {
                // Find the step that produced this cluster.
                // The cluster was the result of a previous merge.
                // Search earlier steps for one whose left∪right = step.left.
                find_producing_step(topo, &step.left, k)
            };
            let right_child = if step.right.len() == 1 {
                nseq + step.right[0]
            } else {
                find_producing_step(topo, &step.right, k)
            };

            // Link parent (node k) ↔ left child
            nodes[k].children[0] = left_child as i32;
            nodes[k].length[0] = step.left_length;
            // Link child back to parent
            let free_slot = find_free_slot(&nodes[left_child]);
            nodes[left_child].children[free_slot] = k as i32;
            nodes[left_child].length[free_slot] = step.left_length;

            // Link parent ↔ right child
            nodes[k].children[1] = right_child as i32;
            nodes[k].length[1] = step.right_length;
            let free_slot = find_free_slot(&nodes[right_child]);
            nodes[right_child].children[free_slot] = k as i32;
            nodes[right_child].length[free_slot] = step.right_length;
        }

        // Compute complement members for internal nodes (direction 2 = "rest").
        let all_seqs: Vec<i32> = (0..nseq as i32).chain(std::iter::once(-1)).collect();
        for k in 0..nseq - 1 {
            let left_set: Vec<i32> = nodes[k].members[0].iter().copied().filter(|&x| x >= 0).collect();
            let right_set: Vec<i32> = nodes[k].members[1].iter().copied().filter(|&x| x >= 0).collect();
            let complement: Vec<i32> = (0..nseq as i32)
                .filter(|s| !left_set.contains(s) && !right_set.contains(s))
                .chain(std::iter::once(-1))
                .collect();
            nodes[k].members[2] = complement;
        }

        // Pre-compute branch weights (calcBranchWeight).
        // For each internal node k and each direction d, compute the branch weight.
        for k in 0..nseq - 1 {
            for d in 0..3 {
                let child_idx = nodes[k].children[d];
                if child_idx >= 0 {
                    let w = calc_w(&nodes, k, child_idx as usize, nseq);
                    nodes[k].branch_weight[d] = w;
                }
            }
        }

        Self { nodes, nseq }
    }

    /// Compute per-sequence weights for a specific branch split.
    /// `step` is the topology step index, `side` is 0 (left) or 1 (right).
    /// Returns a weight vector of length nseq.
    ///
    /// Ports C's `weightFromABranch`.
    pub fn weights_for_branch(&self, topo: &Topology, step: usize, side: usize) -> Vec<f64> {
        let nseq = self.nseq;
        if nseq <= 2 || self.nodes.is_empty() {
            return vec![1.0; nseq];
        }

        let mut result = vec![1.0f64; nseq];

        // Find the split point: btmNode is the subtree for this branch,
        // topNode is the parent.
        let top_node = step;

        // Find which child direction of step matches the requested side.
        let members_to_match = if side == 0 {
            &topo.steps[step].left
        } else {
            &topo.steps[step].right
        };

        let mut btm_dir = 0;
        for d in 0..3 {
            let child = self.nodes[step].children[d];
            if child >= 0 {
                let child_members: Vec<i32> = self.nodes[step].members[d]
                    .iter().copied().filter(|&x| x >= 0).collect();
                let match_members: Vec<i32> = members_to_match.iter().map(|&s| s as i32).collect();
                if child_members == match_members {
                    btm_dir = d;
                    break;
                }
            }
        }

        let btm_node = self.nodes[step].children[btm_dir];
        if btm_node < 0 { return result; }

        // Apply weights recursively from both sides of the split.
        self.weight_from_branch_rec(&mut result, btm_node as usize, top_node);
        self.weight_from_branch_rec(&mut result, top_node, btm_node as usize);

        result
    }

    /// Recursive weight application.
    /// Multiplies branch weights for all sequences reachable through
    /// the children of `node` (excluding the direction toward `from`).
    fn weight_from_branch_rec(&self, result: &mut [f64], node: usize, from: usize) {
        // Leaf: nothing to do (already 1.0).
        if self.is_leaf(node) {
            return;
        }

        // Find directions: which children are NOT the `from` node.
        for d in 0..3 {
            let child = self.nodes[node].children[d];
            if child >= 0 && child as usize != from {
                // Multiply weight for all sequences in this direction.
                let w = self.nodes[node].branch_weight[d];
                for &s in &self.nodes[node].members[d] {
                    if s >= 0 {
                        result[s as usize] *= w;
                    }
                }
                self.weight_from_branch_rec(result, child as usize, node);
            }
        }
    }

    fn is_leaf(&self, idx: usize) -> bool {
        idx >= self.nseq // leaf nodes are at indices nseq..2*nseq-1
    }
}

/// Find which topology step produced a given cluster.
fn find_producing_step(topo: &Topology, cluster: &[usize], before_step: usize) -> usize {
    for k in (0..before_step).rev() {
        let mut merged: Vec<usize> = topo.steps[k].left.clone();
        merged.extend(&topo.steps[k].right);
        merged.sort();
        let mut target = cluster.to_vec();
        target.sort();
        if merged == target {
            return k;
        }
    }
    0 // fallback
}

/// Find a free slot (children[i] == -1) in a node.
fn find_free_slot(node: &WNode) -> usize {
    for i in 0..3 {
        if node.children[i] < 0 { return i; }
    }
    2 // fallback: overwrite last
}

/// Compute the branch weight for a specific direction.
/// Ports C's `calcW`.
fn calc_w(nodes: &[WNode], node_idx: usize, child_idx: usize, nseq: usize) -> f64 {
    const MAX_BW: f64 = 1.0;
    const MIN_BW: f64 = 0.01;

    // Find the direction from node toward child.
    let dir_to_child = (0..3)
        .find(|&d| nodes[node_idx].children[d] == child_idx as i32)
        .unwrap_or(0);

    let top_w = calc_w_single(nodes, node_idx, child_idx, nseq);
    let btm_w = calc_w_single(nodes, child_idx, node_idx, nseq);
    top_w * btm_w
}

/// Compute weight from one side (C's calcW for one node).
fn calc_w_single(nodes: &[WNode], ob_idx: usize, op_idx: usize, nseq: usize) -> f64 {
    const MAX_BW: f64 = 1.0;
    const MIN_BW: f64 = 0.01;

    // Leaf: return 1.0
    if ob_idx >= nseq {
        return 1.0;
    }

    // Find child directions (not toward op).
    let mut dir_ch = Vec::new();
    let mut dir_pa = 0;
    for d in 0..3 {
        if nodes[ob_idx].children[d] == op_idx as i32 {
            dir_pa = d;
        } else if nodes[ob_idx].children[d] >= 0 {
            dir_ch.push(d);
        }
    }
    if dir_ch.len() < 2 {
        return 1.0;
    }

    let a = synthetic_length(nodes, nodes[ob_idx].children[dir_ch[0]] as usize, ob_idx, nseq);
    let b = synthetic_length(nodes, nodes[ob_idx].children[dir_ch[1]] as usize, ob_idx, nseq);
    let c = synthetic_length(nodes, nodes[ob_idx].children[dir_pa] as usize, ob_idx, nseq);

    if c == 0.0 { return MAX_BW; }
    if a == 0.0 || b == 0.0 { return MIN_BW; }

    let s = b * c + c * a + a * b;
    if s == 0.0 { return MAX_BW; }

    let value = a * b * (c + a) * (c + b) / (c * (a + b) * s);
    value.sqrt()
}

/// Compute synthetic length of a subtree.
/// Ports C's `syntheticLength`.
fn synthetic_length(nodes: &[WNode], ob_idx: usize, op_idx: usize, nseq: usize) -> f64 {
    // Find the branch length from ob toward op.
    let len_to_parent = (0..3)
        .find(|&d| nodes[ob_idx].children[d] == op_idx as i32)
        .map(|d| nodes[ob_idx].length[d])
        .unwrap_or(0.0);

    // Leaf: return the branch length.
    if ob_idx >= nseq {
        return len_to_parent;
    }

    // Internal: harmonic mean of children's synthetic lengths, plus branch to parent.
    let mut child_lengths = Vec::new();
    for d in 0..3 {
        let child = nodes[ob_idx].children[d];
        if child >= 0 && child as usize != op_idx {
            child_lengths.push(synthetic_length(nodes, child as usize, ob_idx, nseq));
        }
    }

    if child_lengths.len() < 2 {
        return len_to_parent;
    }

    let (a, b) = (child_lengths[0], child_lengths[1]);
    let value = if a == 0.0 || b == 0.0 {
        0.0
    } else {
        1.0 / (1.0 / a + 1.0 / b) // harmonic mean
    };

    value + len_to_parent
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
