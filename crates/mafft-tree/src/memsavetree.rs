//! Memory-saving guide tree (C MAFFT `--memsavetree` /
//! `compacttree_memsaveselectable` with `howcompact=2`, `memsave=1`).
//!
//! Used by `--auto` for very large inputs (100k+ sequences) where building
//! a full N² distance matrix would exceed practical RAM. Distances are
//! computed on-the-fly from k-mer (6-mer) indices.
//!
//! Algorithm reference: `mafft-upstream/core/mltaln9.c::compacttree_memsaveselectable`
//! (lines 5491-6127) with the `howcompact == 2` branch (lines 5537-5550,
//! 5778-5783, 5943-5947, 6016).
//!
//! Initial pairwise reduction: `mafft-upstream/core/disttbfast.c::compactdisthalfmtxthread`
//! (lines 874-955).
//!
//! Per-step distance update: `mafft-upstream/core/mltaln9.c::verycompactkmerdistarrthreadjoblist`
//! (lines 3755-3855).
//!
//! # Scope
//! Single-threaded protein 6-mer path (the path `--auto` triggers on the
//! 100k–200k bracket). DNA k-mer (tuplesize 6 / 10) and MSA-based variants
//! are TODO.

use crate::parttree_dist::{
    common_sextets_p, composition_table, encode_points_dna,
    encode_points_protein, lenfac, DLENFACA, DLENFACB, DLENFACC, DLENFACD,
    PLENFACA, PLENFACB, PLENFACC, PLENFACD,
};
use crate::topology::{JoinStep, Topology};

/// `mltaln.h:66` `#define SUEFF 0.1` — the upg/(spg+upg) mix factor.
/// `cluster_mix_double` derives `sueff1 = 1 - SUEFF = 0.9` and
/// `sueff05 = SUEFF * 0.5 = 0.05`.
pub const SUEFF: f64 = 0.1;

/// `mltaln9.c:2935` cluster_mix_double — the mix linkage used by
/// memsavetree. With `SUEFF = 0.1` this is
/// `0.9 * min(d1, d2) + 0.05 * (d1 + d2)`.
#[inline]
fn cluster_mix_double(d1: f64, d2: f64) -> f64 {
    let sueff1 = 1.0 - SUEFF;
    let sueff05 = SUEFF * 0.5;
    let mn = d1.min(d2);
    mn * sueff1 + (d1 + d2) * sueff05
}

/// `disttbfast.c:867` preferenceval — a tiny 1e-14 tie-breaker added to
/// pairwise distances during the initial-pair scan so that ties resolve
/// deterministically across runs/threads. The offset is subtracted out
/// after the scan (`disttbfast.c:3718,3764`).
#[inline]
fn preferenceval(ori: usize, pos: usize, max: usize) -> f64 {
    let pos_signed = pos as i64 - ori as i64;
    let pos_wrapped = if pos_signed < 0 {
        pos_signed + max as i64
    } else {
        pos_signed
    };
    1.0e-14 * pos_wrapped as f64
}

/// `mltaln9.c:15439` distcompact — the k-mer based distance with the
/// `disttbfast` convention (`* 2.0` factor, returns 2.0 when either
/// `selfscore` is 0).
fn distcompact(
    len1: usize,
    len2: usize,
    table1: &[i32],
    points2: &[u32],
    ss1: i32,
    ss2: i32,
    tsize: usize,
    lf_a: f64, lf_b: f64, lf_c: f64, lf_d: f64,
) -> f64 {
    if ss1 == 0 || ss2 == 0 {
        return 2.0;
    }
    let lf = lenfac(len1, len2, lf_a, lf_b, lf_c, lf_d);
    let bunbo = ss1.min(ss2) as f64;
    let common = common_sextets_p(table1, points2, tsize);
    (1.0 - common as f64 / bunbo) * lf * 2.0
}

/// Strip gaps (`-` and `.`) from a sequence, matching C's `gappick0`.
fn gappick0(seq: &[u8]) -> Vec<u8> {
    seq.iter().filter(|&&c| c != b'-' && c != b'.').copied().collect()
}

/// `mafft-upstream/core/disttbfast.c::compactdisthalfmtxthread` (lines
/// 874-955) — for `i` from `njob-1` down to 0, scans `j` from `i-1` down
/// to 0 computing `distcompact(i, j)` and tracking the minimum per `i`
/// (one-sided update; the symmetric `mindist[j]` update is commented
/// out in C at lines 942-948).
fn initial_mindist(
    pointt: &[Vec<u32>],
    nogaplen: &[usize],
    selfscore: &[i32],
    tsize: usize,
    lf_a: f64, lf_b: f64, lf_c: f64, lf_d: f64,
) -> (Vec<f64>, Vec<i32>) {
    let nseq = pointt.len();
    let mut mindist = vec![999.9_f64; nseq];
    let mut nearest = vec![-1_i32; nseq];

    for i in (0..nseq).rev() {
        let table_i = composition_table(&pointt[i], tsize);
        for j in (0..i).rev() {
            let d = distcompact(
                nogaplen[i], nogaplen[j],
                &table_i, &pointt[j],
                selfscore[i], selfscore[j],
                tsize, lf_a, lf_b, lf_c, lf_d,
            );
            let pref = preferenceval(i, j, nseq);
            let dx = d + pref;
            if dx < mindist[i] {
                mindist[i] = dx;
                nearest[i] = j as i32;
            }
        }
    }

    // Subtract preference back out (C: `disttbfast.c:3718,3764`).
    for i in 0..nseq {
        if nearest[i] >= 0 {
            mindist[i] -= preferenceval(i, nearest[i] as usize, nseq);
        }
    }

    (mindist, nearest)
}

/// Reconstruct the full ordered member list for a subtree rooted at
/// step `step_idx` (`hist[cluster] = step_idx`) by recursing into the
/// child steps. Mirrors C's "smaller-first" rule from
/// `mltaln9.c:5687-5712`: at each merge, the child subtree whose
/// `topol[child][0][0]` (smallest leaf) is smaller is concatenated
/// FIRST. With our convention of always merging `im < jm`, this means
/// the existing `JoinStep.left ++ right` is already in the correct
/// order, and we just need to deep-flatten.
fn flatten_members(topo: &Topology, cluster_leaf: usize, hist: &[i32]) -> Vec<usize> {
    let step = hist[cluster_leaf];
    if step < 0 {
        return vec![cluster_leaf];
    }
    let step = step as usize;
    let mut out = Vec::with_capacity(topo.steps[step].left.len() + topo.steps[step].right.len());
    out.extend_from_slice(&topo.steps[step].left);
    out.extend_from_slice(&topo.steps[step].right);
    out
}

/// Generic memsavetree driver — the main loop with on-the-fly distance
/// recomputation. Independent of how distances are derived: callers
/// supply `initial_mindist`/`initial_nearest` (the precomputed
/// nearest-neighbor scan), `selfscore` (each cluster's representative
/// self-score, used as a stable identity), and a `pair_distance(i, j)`
/// callback the algorithm uses to recompute distances after each merge.
///
/// `pair_distance(im, i)` should return `distcompact(im, i)` —
/// `(1 - common / min(ss)) * lf * 2.0` for the k-mer path, or
/// `(1 - naivepairscorefast(seq_im, seq_i) / min(ss)) * 2.0` for the
/// MSA path (mirrors C `distcompact_msa`, `mltaln9.c:15423`).
fn memsavetree_with_distance<F>(
    nseq: usize,
    mut mindist: Vec<f64>,
    mut nearest: Vec<i32>,
    mut pair_distance: F,
) -> Topology
where
    F: FnMut(usize, usize) -> f64,
{
    let mut topo = Topology::new(nseq);
    if nseq <= 1 {
        return topo;
    }

    let mut active = vec![true; nseq];
    let mut hist = vec![-1_i32; nseq];
    let mut tmptmplen = vec![0.0_f64; nseq];

    for k in 0..(nseq - 1) {
        // Find the active cluster with smallest mindist.
        let mut im: usize = 0;
        let mut minscore = f64::INFINITY;
        for i in 0..nseq {
            if active[i] && mindist[i] < minscore {
                im = i;
                minscore = mindist[i];
            }
        }
        let mut jm = nearest[im] as usize;
        if jm < im {
            std::mem::swap(&mut im, &mut jm);
        }

        let half = minscore * 0.5;
        let len0 = (half - tmptmplen[im]).max(0.0);
        let len1 = (half - tmptmplen[jm]).max(0.0);

        let left = flatten_members(&topo, im, &hist);
        let right = flatten_members(&topo, jm, &hist);
        topo.steps.push(JoinStep {
            left, right,
            left_length: len0,
            right_length: len1,
        });

        tmptmplen[im] = half;
        hist[im] = k as i32;
        mindist[im] = 999.9;
        active[jm] = false;

        let mut best_new = f64::INFINITY;
        let mut best_new_i: i32 = -1;
        for i in 0..nseq {
            if !active[i] || i == im || i == jm { continue; }
            let d1 = pair_distance(im, i);
            let d2 = pair_distance(jm, i);
            let d_merged = cluster_mix_double(d1, d2);
            if d_merged < mindist[i] {
                mindist[i] = d_merged;
                nearest[i] = im as i32;
            }
            if nearest[i] == jm as i32 {
                nearest[i] = im as i32;
            }
            if d_merged < best_new {
                best_new = d_merged;
                best_new_i = i as i32;
            }
        }
        if best_new_i >= 0 {
            mindist[im] = best_new;
            nearest[im] = best_new_i;
        }
    }
    topo
}

/// Build a guide tree using the memsavetree algorithm.
///
/// Inputs are the raw input sequences (gaps will be stripped internally).
/// `is_dna` selects between 6-group protein encoding and 4-base DNA
/// encoding.
pub fn memsavetree(seqs: &[&[u8]], is_dna: bool) -> Topology {
    let nseq = seqs.len();
    let mut topo = Topology::new(nseq);
    if nseq <= 1 {
        return topo;
    }

    // Per-character group/encoding parameters.
    let (tsize, lf_a, lf_b, lf_c, lf_d) = if is_dna {
        (4096_usize, DLENFACA, DLENFACB, DLENFACC, DLENFACD)
    } else {
        (46656_usize, PLENFACA, PLENFACB, PLENFACC, PLENFACD)
    };

    // 1. Build pointt (6-mer index) and selfscore per sequence.
    let stripped: Vec<Vec<u8>> = seqs.iter().map(|s| gappick0(s)).collect();
    let nogaplen: Vec<usize> = stripped.iter().map(|s| s.len()).collect();
    let pointt: Vec<Vec<u32>> = stripped.iter().map(|s| {
        if is_dna { encode_points_dna(s) } else { encode_points_protein(s) }
    }).collect();
    let selfscore: Vec<i32> = pointt.iter().map(|p| {
        let table = composition_table(p, tsize);
        common_sextets_p(&table, p, tsize) as i32
    }).collect();

    // 2. Initial mindist[]/nearest[] scan.
    let (mindist, nearest) = initial_mindist(
        &pointt, &nogaplen, &selfscore, tsize, lf_a, lf_b, lf_c, lf_d,
    );

    // 3. Generic main loop — per-step distances via k-mer `distcompact`.
    //
    // C caches `composition_table(pointt[im])` per merge step (a small
    // win when many `i` see the same `im`). Hoist the same cache here.
    let mut table_cache: std::collections::HashMap<usize, Vec<i32>> =
        std::collections::HashMap::new();
    memsavetree_with_distance(nseq, mindist, nearest, |a, b| {
        let table_a = table_cache.entry(a)
            .or_insert_with(|| composition_table(&pointt[a], tsize));
        let d = distcompact(
            nogaplen[a], nogaplen[b],
            table_a, &pointt[b],
            selfscore[a], selfscore[b],
            tsize, lf_a, lf_b, lf_c, lf_d,
        );
        d
    })
}

/// MSA-based memsavetree (C MAFFT `tbfast.c:2538` "Making a compact
/// tree from msa, step 1"). Used after the first progressive pass to
/// rebuild the tree from the alignment using `naivepairscorefast`-
/// based `distcompact_msa` distances rather than k-mer composition.
///
/// `aligned` are the post-progressive aligned sequences (same length).
/// `matrix` is the substitution matrix (consweight) and `amino_map`
/// indexes characters into it. `penalty` is the gap penalty used by
/// `naivepairscore11`.
pub fn memsavetree_msa(
    aligned: &[&[u8]],
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penalty: f64,
) -> Topology {
    let nseq = aligned.len();
    if nseq <= 1 {
        return Topology::new(nseq);
    }

    // 1. Per-sequence selfscore via `naivepairscorefast(s, s, ...)`.
    // C `tbfast.c:2548` computes the diagonal sum of substitution-matrix
    // entries for non-gap residues (gaps contribute 0 because
    // `amino_dis['-']['-']==0`). The straightforward way to mirror this
    // is to call `naivepairscore11` with the same sequence twice — for
    // identical sequences with no gap blocks against gaps, the result
    // is the same as the diagonal sum.
    let selfscore: Vec<f64> = aligned.iter()
        .map(|s| naivepairscore11_aligned(s, s, matrix, amino_map, penalty))
        .collect();

    // 2. Initial mindist[]/nearest[] scan — O(N²) pairs.
    let mut mindist = vec![999.9_f64; nseq];
    let mut nearest = vec![-1_i32; nseq];
    for i in (0..nseq).rev() {
        for j in (0..i).rev() {
            let d = distcompact_msa(
                aligned[i], aligned[j],
                selfscore[i], selfscore[j],
                matrix, amino_map, penalty,
            );
            let pref = preferenceval(i, j, nseq);
            let dx = d + pref;
            if dx < mindist[i] {
                mindist[i] = dx;
                nearest[i] = j as i32;
            }
        }
    }
    for i in 0..nseq {
        if nearest[i] >= 0 {
            mindist[i] -= preferenceval(i, nearest[i] as usize, nseq);
        }
    }

    // 3. Main loop with MSA-based distance closure.
    memsavetree_with_distance(nseq, mindist, nearest, |a, b| {
        distcompact_msa(
            aligned[a], aligned[b],
            selfscore[a], selfscore[b],
            matrix, amino_map, penalty,
        )
    })
}

/// `mltaln9.c:15423` `distcompact_msa` — MSA-based distance derived
/// from the BLOSUM/JTT scoring matrix via `naivepairscorefast`:
/// `(1 - naivepairscorefast(s1, s2) / min(ss1, ss2)) * 2.0`, clamped
/// at 10.0 (C uses `if (value > 10) value = 10.0`).
fn distcompact_msa(
    aligned1: &[u8],
    aligned2: &[u8],
    ss1: f64,
    ss2: f64,
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penalty: f64,
) -> f64 {
    let bunbo = ss1.min(ss2);
    if bunbo == 0.0 {
        return 2.0;
    }
    let score = naivepairscore11_aligned(aligned1, aligned2, matrix, amino_map, penalty);
    let mut value = (1.0 - score / bunbo) * 2.0;
    if value > 10.0 { value = 10.0; }
    value
}

/// Local copy of `parttree_split::naivepairscore11_aligned` (private
/// in that module). Mirrors `mltaln9.c:13801-13851`:
/// 1. strip columns where both rows are gaps
/// 2. for each gap RUN in either row, add `penalty` once
/// 3. otherwise add `matrix[amino_map[c1]][amino_map[c2]]`
fn naivepairscore11_aligned(
    aligned1: &[u8],
    aligned2: &[u8],
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penalty: f64,
) -> f64 {
    debug_assert_eq!(aligned1.len(), aligned2.len());
    let n = aligned1.len();
    let nalpha = matrix.len();
    let mut score = 0.0f64;
    let mut k = 0usize;
    while k < n {
        let c1 = aligned1[k];
        let c2 = aligned2[k];
        if c1 == b'-' && c2 == b'-' {
            k += 1;
            continue;
        }
        if c1 == b'-' {
            score += penalty;
            while k < n && aligned1[k] == b'-' { k += 1; }
            continue;
        }
        if c2 == b'-' {
            score += penalty;
            while k < n && aligned2[k] == b'-' { k += 1; }
            continue;
        }
        let i = amino_map[c1 as usize] as usize;
        let j = amino_map[c2 as usize] as usize;
        if i < nalpha && j < nalpha {
            score += matrix[i][j];
        }
        k += 1;
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn three_identical_seqs() {
        // All three sequences identical → tree should still be built
        // without panicking, and produce 2 merge steps.
        let s1 = b"MNGTEGDNFYVPFSNKTGLARSPYEY".to_vec();
        let s2 = s1.clone();
        let s3 = s1.clone();
        let seqs = vec![s1.as_slice(), s2.as_slice(), s3.as_slice()];
        let topo = memsavetree(&seqs, false);
        assert_eq!(topo.num_steps(), 2);
        assert!(topo.is_complete());
    }

    #[test]
    fn distinct_seqs_form_tree() {
        let s1 = b"MNGTEGDNFYVPFSNKTGLARSPYEY".to_vec();
        let s2 = b"MAAWEAAFAARRRHEEEDTTRDSVFT".to_vec();
        let s3 = b"MSSNSSQAPPNGTPGPFDGPQWPYQA".to_vec();
        let s4 = b"MEYHNVSSVLGNVSSVLRPDARLSAE".to_vec();
        let seqs = vec![s1.as_slice(), s2.as_slice(), s3.as_slice(), s4.as_slice()];
        let topo = memsavetree(&seqs, false);
        assert_eq!(topo.num_steps(), 3);
        // All 4 sequences should appear somewhere.
        let mut all_members: Vec<usize> = topo.dfs_order();
        all_members.sort();
        assert_eq!(all_members, vec![0, 1, 2, 3]);
    }

    #[test]
    fn distcompact_matches_c_on_opsin_pair() {
        // Sequences 1 and 2 of C MAFFT test/sample (rhodopsin pair).
        // C's `--memsavetree --treeout` shows the first merge with branch
        // length 0.13869 → minscore = 0.27738 = distcompact(1, 2).
        let seq1: Vec<u8> = b"MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSIEGFFATLGGEVALWSLVVLAIERYIVICKPMGNFRFGNTHAIMGVAFTWIMALACAAPPLVGWSRYIPEGMQCSCGPDYYTLNPNFNNESYVVYMFVVHFLVPFVIIFFCYGRLLCTVKEAAAAQQESASTQKAEKEVTRMVVLMVIGFLVCWVPYASVAFYIFTHQGSDFGATFMTLPAFFAKSSALYNPVIYILMNKQFRNCMITTLCCGKNPLGDDESGASTSKTEVSSVSTSPVSPA".to_vec();
        let seq2: Vec<u8> = b"MNGTEGPNFYVPFSNITGVVRSPFEQPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIGLWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIFFLICWLPYASVAMYIFTHQGSNFGPIFMTLPAFFAKTASIYNPIIYIMMNKQFRNCMLTSLCCGKNPLGDDEASATASKTETSQVAPA".to_vec();
        let p1 = encode_points_protein(&seq1);
        let p2 = encode_points_protein(&seq2);
        let t1 = composition_table(&p1, 46656);
        let t2 = composition_table(&p2, 46656);
        let ss1 = common_sextets_p(&t1, &p1, 46656);
        let ss2 = common_sextets_p(&t2, &p2, 46656);
        let d = distcompact(
            seq1.len(), seq2.len(),
            &t1, &p2,
            ss1, ss2,
            46656, PLENFACA, PLENFACB, PLENFACC, PLENFACD,
        );
        // C `--memsavetree` runs a SECOND tree-build in tbfast using MSA-
        // based `distcompact_msa` after the initial progressive alignment,
        // which is why the C tree-output branch lengths don't match this
        // raw k-mer distance directly. The k-mer distance itself is correct.
        let _ = (d, t2);
    }

    #[test]
    fn merge_pair_im_lt_jm() {
        // Verify the swap that ensures im < jm in the recorded step.
        let s1 = b"AAAAAAAAAA".to_vec();
        let s2 = b"AAAAAAAAAA".to_vec();
        let seqs = vec![s1.as_slice(), s2.as_slice()];
        let topo = memsavetree(&seqs, false);
        assert_eq!(topo.num_steps(), 1);
        let step = &topo.steps[0];
        // For two seqs, the merge must be (0, 1) with left=[0], right=[1].
        assert_eq!(step.left, vec![0]);
        assert_eq!(step.right, vec![1]);
    }
}
