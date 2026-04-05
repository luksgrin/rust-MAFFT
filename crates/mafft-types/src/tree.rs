/// A node in a bifurcating guide tree.
///
/// Uses index-based references (into a `Vec<TreeNode>`) instead of raw
/// pointers. `None` means this is a leaf.
#[derive(Debug, Clone)]
pub struct TreeNode {
    /// Left child index, or `None` for a leaf.
    pub left: Option<usize>,
    /// Right child index, or `None` for a leaf.
    pub right: Option<usize>,
    /// Branch length to left child.
    pub left_length: f64,
    /// Branch length to right child.
    pub right_length: f64,
    /// Sequence indices that belong to this subtree.
    pub members: Vec<usize>,
}

/// A guide tree built by neighbor-joining or UPGMA.
///
/// Nodes are stored in a flat arena. The root is the last node.
/// Replaces the C `int ***topol` + `double **len` representation.
#[derive(Debug, Clone)]
pub struct GuideTree {
    pub nodes: Vec<TreeNode>,
}

impl GuideTree {
    /// Create an empty tree with pre-allocated capacity.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(2 * n),
        }
    }

    /// The root node index (last node in the arena).
    pub fn root(&self) -> Option<usize> {
        if self.nodes.is_empty() {
            None
        } else {
            Some(self.nodes.len() - 1)
        }
    }

    /// Number of leaves (sequences) in the tree.
    pub fn num_leaves(&self) -> usize {
        self.nodes.iter().filter(|n| n.left.is_none() && n.right.is_none()).count()
    }
}

/// Tree dependency information for memory-efficient computation.
///
/// Replaces C `Treedep` struct. Used during progressive alignment to
/// track which nodes have been processed.
#[derive(Debug, Clone)]
pub struct TreeDependency {
    pub child0: i32,
    pub child1: i32,
    pub done: bool,
    pub dist_from_tip: f64,
}

impl Default for TreeDependency {
    fn default() -> Self {
        Self {
            child0: -1,
            child1: -1,
            done: false,
            dist_from_tip: 0.0,
        }
    }
}
