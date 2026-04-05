/// Production UPGMA tree builder with nearest-neighbor heuristic.
///
/// Ports the C `fixed_musclesupg_double_realloc_nobk_halfmtx()` from mltaln9.c.
///
/// This is the tree builder used in practice by MAFFT. It uses a
/// nearest-neighbor cache for O(n²) performance, a doubly-linked list
/// for efficient cluster removal, and configurable distance update functions.

use crate::distance::DistanceMatrix;
use crate::topology::{JoinStep, Topology};

/// Cluster distance update strategy.
#[derive(Debug, Clone, Copy)]
pub enum ClusterMethod {
    /// Average linkage (UPGMA): `(d1 + d2) / 2`
    Average,
    /// Single linkage: `min(d1, d2)`
    Minimum,
    /// Mixed: weighted combination of min and average.
    /// `min(d1,d2) * (1-sueff) + (d1+d2)/2 * sueff`
    Mix { sueff: f64 },
}

impl ClusterMethod {
    fn compute(&self, d1: f64, d2: f64) -> f64 {
        match self {
            Self::Average => (d1 + d2) * 0.5,
            Self::Minimum => d1.min(d2),
            Self::Mix { sueff } => {
                d1.min(d2) * (1.0 - sueff) + (d1 + d2) * 0.5 * sueff
            }
        }
    }
}

impl Default for ClusterMethod {
    fn default() -> Self {
        Self::Mix { sueff: 0.1 } // MAFFT default: SUEFF = 0.1
    }
}

/// Build a guide tree using the production MUSCLE-style UPGMA with
/// nearest-neighbor heuristic.
///
/// This is significantly faster than textbook NJ for large inputs
/// because it caches each cluster's nearest neighbor and only recomputes
/// when necessary.
pub fn musclesupg(dist: &DistanceMatrix, method: ClusterMethod) -> Topology {
    let n = dist.nseq;
    if n <= 1 {
        return Topology::new(n);
    }
    if n == 2 {
        let d = dist.get(0, 1);
        let mut topo = Topology::new(2);
        topo.steps.push(JoinStep {
            left: vec![0],
            right: vec![1],
            left_length: d * 0.5,
            right_length: d * 0.5,
        });
        return topo;
    }

    // Working copy: half-matrix stored as Vec<Vec<f64>>
    // eff[i] has entries for j > i: eff[i][j-i-1] (same as DistanceMatrix)
    let mut eff: Vec<Option<Vec<f64>>> = (0..n)
        .map(|i| {
            let row: Vec<f64> = (0..n - i - 1).map(|j| dist.get(i, i + j + 1)).collect();
            Some(row)
        })
        .collect();

    // Active cluster chain (simulates Bchain doubly-linked list)
    let mut active: Vec<bool> = vec![true; n];
    let mut next: Vec<Option<usize>> = (0..n).map(|i| if i + 1 < n { Some(i + 1) } else { None }).collect();
    let mut prev: Vec<Option<usize>> = (0..n).map(|i| if i > 0 { Some(i - 1) } else { None }).collect();
    let first_active = 0usize;

    // Nearest-neighbor cache
    let mut mindisfrom = vec![f64::MAX; n];
    let mut nearest = vec![0usize; n];

    // Initialize nearest neighbors
    for i in 0..n {
        find_nearest(i, &eff, &active, n, &mut mindisfrom[i], &mut nearest[i]);
    }

    // Cluster tracking
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut tmplen = vec![0.0f64; n];

    let mut topo = Topology::new(n);

    for _step in 0..(n - 1) {
        // Find cluster with minimum distance to its nearest neighbor
        let mut min_score = f64::MAX;
        let mut im = 0;

        let mut idx = Some(first_active);
        while let Some(i) = idx {
            if active[i] && mindisfrom[i] < min_score {
                min_score = mindisfrom[i];
                im = i;
            }
            idx = next[i];
        }

        let jm = nearest[im];
        if !active[jm] {
            // Nearest was invalidated, recompute
            find_nearest(im, &eff, &active, n, &mut mindisfrom[im], &mut nearest[im]);
            continue; // retry with updated nearest
        }

        // Ensure im < jm for consistent half-matrix access
        let (im, jm) = if im < jm { (im, jm) } else { (jm, im) };

        // Branch lengths (UPGMA-style)
        let node_height = min_score * 0.5;
        let left_len = node_height - tmplen[im];
        let right_len = node_height - tmplen[jm];

        topo.steps.push(JoinStep {
            left: members[im].clone(),
            right: members[jm].clone(),
            left_length: left_len,
            right_length: right_len,
        });

        // Update distances: merge jm into im
        mindisfrom[im] = f64::MAX;

        let mut idx = Some(first_active);
        while let Some(i) = idx {
            if !active[i] || i == im || i == jm {
                idx = next[i];
                continue;
            }

            let d_im = get_half(&eff, i, im);
            let d_jm = get_half(&eff, i, jm);
            let new_dist = method.compute(d_im, d_jm);

            set_half(&mut eff, i, im, new_dist);

            // Update nearest cache for i
            if new_dist < mindisfrom[i] {
                mindisfrom[i] = new_dist;
                nearest[i] = im;
            } else if nearest[i] == jm {
                nearest[i] = im;
                if d_jm < d_im {
                    // jm was closer, recompute
                    find_nearest(i, &eff, &active, n, &mut mindisfrom[i], &mut nearest[i]);
                }
            }

            // Update nearest for im
            if new_dist < mindisfrom[im] {
                mindisfrom[im] = new_dist;
                nearest[im] = i;
            }

            idx = next[i];
        }

        // Merge members
        let jm_members = std::mem::take(&mut members[jm]);
        members[im].extend(jm_members);
        tmplen[im] = node_height;

        // Remove jm from chain
        active[jm] = false;
        if let Some(p) = prev[jm] {
            next[p] = next[jm];
        }
        if let Some(nx) = next[jm] {
            prev[nx] = prev[jm];
        }
        // Free jm's distance row
        eff[jm] = None;
    }

    topo
}

/// Get distance from half-matrix, handling i < j vs i > j.
fn get_half(eff: &[Option<Vec<f64>>], i: usize, j: usize) -> f64 {
    if i == j { return 0.0; }
    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
    eff[lo].as_ref().map_or(f64::MAX, |row| {
        row.get(hi - lo - 1).copied().unwrap_or(f64::MAX)
    })
}

/// Set distance in half-matrix.
fn set_half(eff: &mut [Option<Vec<f64>>], i: usize, j: usize, val: f64) {
    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
    if let Some(row) = &mut eff[lo] {
        if hi - lo - 1 < row.len() {
            row[hi - lo - 1] = val;
        }
    }
}

/// Find nearest active neighbor for cluster `pos`.
fn find_nearest(
    pos: usize,
    eff: &[Option<Vec<f64>>],
    active: &[bool],
    n: usize,
    out_mindis: &mut f64,
    out_nearest: &mut usize,
) {
    let mut best_dist = f64::MAX;
    let mut best_idx = 0;

    for j in 0..n {
        if j == pos || !active[j] { continue; }
        let d = get_half(eff, pos, j);
        if d < best_dist {
            best_dist = d;
            best_idx = j;
        }
    }

    *out_mindis = best_dist;
    *out_nearest = best_idx;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn musclesupg_two_sequences() {
        let mut dm = DistanceMatrix::new(2);
        dm.set(0, 1, 0.4);
        let topo = musclesupg(&dm, ClusterMethod::Average);
        assert_eq!(topo.num_steps(), 1);
        assert!((topo.steps[0].left_length - 0.2).abs() < 1e-10);
    }

    #[test]
    fn musclesupg_joins_closest() {
        let mut dm = DistanceMatrix::new(4);
        dm.set(0, 1, 0.1); // closest
        dm.set(0, 2, 0.5);
        dm.set(0, 3, 0.6);
        dm.set(1, 2, 0.5);
        dm.set(1, 3, 0.6);
        dm.set(2, 3, 0.2); // second closest

        let topo = musclesupg(&dm, ClusterMethod::Average);
        assert!(topo.is_complete());

        // First join should be 0+1 (distance 0.1)
        let first = &topo.steps[0];
        let mut joined: Vec<usize> = first.left.iter().chain(first.right.iter()).copied().collect();
        joined.sort();
        assert_eq!(joined, vec![0, 1]);
    }

    #[test]
    fn musclesupg_preserves_all() {
        let mut dm = DistanceMatrix::new(5);
        for i in 0..5 {
            for j in (i + 1)..5 {
                dm.set(i, j, (i + j + 1) as f64 * 0.1);
            }
        }
        let topo = musclesupg(&dm, ClusterMethod::default());
        let last = topo.steps.last().unwrap();
        let mut all: Vec<usize> = last.left.iter().chain(last.right.iter()).copied().collect();
        all.sort();
        assert_eq!(all, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn mix_method() {
        let m = ClusterMethod::Mix { sueff: 0.1 };
        // min(2,3) * 0.9 + (2+3)/2 * 0.1 = 1.8 + 0.25 = 2.05
        let result = m.compute(2.0, 3.0);
        assert!((result - 2.05).abs() < 1e-10);
    }
}
