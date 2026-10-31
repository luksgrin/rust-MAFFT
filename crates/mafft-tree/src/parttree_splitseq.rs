//! Full recursive port of `splittbfast.c::splitseq_mq` (lines 1174-2642).
//!
//! `splitseq_mq` partitions its `nin` sequences around up to `picksize`
//! pivots, drops redundant pivots, assigns every sequence to its closest
//! surviving pivot ("yuko"), recurses into each yuko's member group, and
//! finally merges the yuko groups along a UPGMA tree built on the
//! pivot-to-pivot distances (`pairalign`, `splittbfast.c:2444-2554`).
//! Every group with more than one distinct sequence is therefore aligned
//! internally (by its own recursive call) before the parent merges it.
//!
//! The run produces the three things C derives from one `splittbfast`
//! call:
//! - the merge steps in C execution order (children first, in `treeorder`,
//!   then the level's own UPGMA steps), each oriented as C's
//!   `if( *mem1 < *mem2 ) pairalign( mem1, mem2 ) else pairalign( mem2, mem1 )`
//!   with both member lists sorted ascending (`qsort( mem, intcompare )`);
//! - the leaf-emission `order[]` used for `--reorder`;
//! - the `--treeout` string.
//!
//! Distances are supplied by a [`SplitScorer`], so the same recursion
//! serves the k-tuple pass (`--parttree` CALL 1), the `-Z` / `fromaln`
//! pass (`--parttree` CALL 2+) and the `-U` DP pass (`--dpparttree`).
//!
//! Pivot selection consumes the C library `rand()` (via `rnd()`,
//! `mltaln.h:56`) once `nin > picksize`. `splittbfast` never seeds it, so
//! each call starts from seed 1. [`BsdRand`] reproduces the macOS / FreeBSD
//! `rand()` sequence, the same platform pinning as the BSD `qsort` port
//! used for `dcompare`.

use std::cell::RefCell;

use crate::distance::DistanceMatrix;
use crate::musclesupg::{musclesupg, ClusterMethod};
use crate::parttree_dist::{common_sextets_p, composition_table, MAX6DIST};
use crate::parttree_split::c_normalized_subtree;
use crate::topology::{JoinStep, Topology};

/// `splittbfast.c::rnd()` source: macOS / FreeBSD libc `rand()`
/// (Park-Miller "minimal standard", multiplier 16807, modulus 2^31 - 1),
/// unseeded state = 1.
#[derive(Debug, Clone)]
pub struct BsdRand {
    next: u64,
}

impl Default for BsdRand {
    fn default() -> Self {
        Self { next: 1 }
    }
}

impl BsdRand {
    pub fn rand(&mut self) -> u32 {
        if self.next == 0 {
            self.next = 123_459_876;
        }
        let ctx = self.next as i64;
        let hi = ctx / 127_773;
        let lo = ctx % 127_773;
        let mut x = 16_807 * lo - 2_836 * hi;
        if x < 0 {
            x += 0x7fff_ffff;
        }
        self.next = x as u64;
        (x as u64 % (0x7fff_ffffu64 + 1)) as u32
    }

    /// `rnd()` = `( 1.0 / ( RAND_MAX + 1.0 ) ) * rand()`.
    pub fn rnd(&mut self) -> f64 {
        (1.0 / (2_147_483_647.0 + 1.0)) * self.rand() as f64
    }
}

/// Distance callback for one `splittbfast` mode. `dist(r, o)` is the value
/// C stores when `r` plays the role of the reference / pivot / yuko and `o`
/// the other sequence (both original input indices), including the mode's
/// clamping. `same(a, b)` is C's duplicate test
/// (`shimon` + `strcmp` on the gap-free sequences).
pub trait SplitScorer {
    fn dist(&self, r: usize, o: usize) -> f64;
    fn same(&self, a: usize, b: usize) -> bool;
}

/// Per-sequence `Scores` fields used by `splitseq_mq`.
#[derive(Debug, Clone)]
pub struct SplitEntry {
    pub numinseq: usize,
    pub selfscore: i64,
    pub orilen: usize,
    pub score: f64,
}

/// Configuration shared by every recursion level of one `splittbfast` call.
#[derive(Debug, Clone, Copy)]
pub struct SplitParams {
    /// `-p` (`picksize`, default `PICKSIZE` = 50).
    pub picksize: usize,
    /// Effective `tokyoripara` (`splittbfast.c:2748-2754`: 0.70, or 0.0
    /// when `picksize > njob`).
    pub tokyoripara: f64,
    /// `-Z`: zero `maxdist` when the closest and farthest entries are the
    /// same sequence (`splittbfast.c:1461-1468`).
    pub fromaln: bool,
}

impl SplitParams {
    /// Defaults the mafft script uses (`-p 50`), with C's
    /// `picksize > njob → tokyoripara = 0` rule applied for `njob`.
    pub fn for_njob(njob: usize, fromaln: bool) -> Self {
        let picksize = crate::parttree_dist::PICKSIZE;
        let tokyoripara = if picksize > njob { 0.0 } else { crate::parttree_dist::TOKYORIPARA };
        Self { picksize, tokyoripara, fromaln }
    }
}

/// Result of one `splittbfast` call.
#[derive(Debug, Clone)]
pub struct SplitRun {
    /// Merge steps in C execution order, oriented and sorted as C's
    /// `pairalign` receives them.
    pub topology: Topology,
    /// Leaf-emission order (`*orderpos++ = scores[j].numinseq`).
    pub order: Vec<usize>,
    /// The `--treeout` string as `splittbfast` writes it (without the
    /// trailing newline added by `fprintf( treefp, "%s\n", *tree )`).
    pub tree: String,
}

struct State<'a, S: SplitScorer + ?Sized> {
    scorer: &'a S,
    params: SplitParams,
    rng: BsdRand,
    steps: Vec<JoinStep>,
    order: Vec<usize>,
}

/// Run `splitseq_mq` on `entries` (the `scores[]` array `main` builds, in
/// input order) and collect the merge steps, order and tree.
pub fn run_splitseq<S: SplitScorer + ?Sized>(
    entries: Vec<SplitEntry>,
    scorer: &S,
    params: SplitParams,
) -> SplitRun {
    let nseq = entries.len();
    let mut st = State {
        scorer,
        params,
        rng: BsdRand::default(),
        steps: Vec::with_capacity(nseq.saturating_sub(1)),
        order: Vec::with_capacity(nseq),
    };
    let tree = splitseq_mq(&mut st, entries, false, None);
    let mut topology = Topology::new(nseq);
    topology.steps = st.steps;
    SplitRun { topology, order: st.order, tree }
}

/// C's `topol[l][0]` / `topol[l][1]` member lists for every UPGMA step:
/// each side in the smaller-first-element normalized order of
/// `mltaln9.c:8184-8197`, and side 0 the cluster holding the smaller
/// representative (`im < jm`).
fn c_topol(steps: &[JoinStep]) -> Vec<(Vec<usize>, Vec<usize>)> {
    steps
        .iter()
        .map(|s| {
            let a = c_normalized_subtree(steps, &s.left);
            let b = c_normalized_subtree(steps, &s.right);
            if a[0] < b[0] { (a, b) } else { (b, a) }
        })
        .collect()
}

fn leaf_tree(scores: &[SplitEntry]) -> String {
    let nin = scores.len();
    let mut t = String::from("\n");
    if nin > 1 {
        t.push('(');
    }
    for (j, s) in scores.iter().enumerate() {
        t.push_str(&(s.numinseq + 1).to_string());
        if j + 1 < nin {
            t.push(',');
        }
    }
    t.push_str(if nin > 1 { ")\n" } else { "\n" });
    t
}

fn splitseq_mq<S: SplitScorer + ?Sized>(
    st: &mut State<'_, S>,
    mut scores: Vec<SplitEntry>,
    uniform_leaf: bool,
    qinoya: Option<usize>,
) -> String {
    let nin = scores.len();
    if nin == 0 {
        return String::new();
    }
    if nin < 2 || uniform_leaf {
        st.order.extend(scores.iter().map(|s| s.numinseq));
        return leaf_tree(&scores);
    }

    // `uselongest`: first entry with the maximum selfscore becomes the
    // reference (`splittbfast.c:1316-1343`).
    let mut belongto = 0usize;
    let mut selfscore0 = scores[0].selfscore;
    for (i, s) in scores.iter().enumerate() {
        if s.selfscore > selfscore0 {
            selfscore0 = s.selfscore;
            belongto = i;
        }
    }
    if belongto != 0 {
        scores.swap(0, belongto);
    }

    // Distances from the reference; skipped when it is the parent's
    // reference, whose distances the entries already carry (`:1360`).
    if qinoya != Some(scores[0].numinseq) {
        let r = scores[0].numinseq;
        for s in scores.iter_mut() {
            s.score = st.scorer.dist(r, s.numinseq);
        }
    }

    crate::bsd_qsort::bsd_qsort(&mut scores, |a, b| {
        use std::cmp::Ordering;
        if a.score > b.score { return Ordering::Greater; }
        if a.score < b.score { return Ordering::Less; }
        if a.selfscore < b.selfscore { return Ordering::Greater; }
        if a.selfscore > b.selfscore { return Ordering::Less; }
        if a.orilen < b.orilen { return Ordering::Greater; }
        if a.orilen > b.orilen { return Ordering::Less; }
        Ordering::Equal
    });

    let mut maxdist = scores[nin - 1].score;
    if st.params.fromaln && st.scorer.same(scores[0].numinseq, scores[nin - 1].numinseq) {
        maxdist = 0.0;
    }

    // Pivot selection (`splittbfast.c:1491-1574`).
    let same_at = |st: &State<'_, S>, a: usize, b: usize| {
        st.scorer.same(scores[a].numinseq, scores[b].numinseq)
    };
    let mut pickkouho: Vec<usize> = (1..=nin).collect();
    let mut nkouho = nin - 1;
    let mut picks: Vec<usize> = vec![0];
    if nkouho > 0 {
        let picktmp = pickkouho[nkouho - 1];
        nkouho -= 1;
        if !same_at(st, picktmp, 0) {
            picks.push(picktmp);
        }
    }
    let mut first = true;
    while picks.len() < st.params.picksize && nkouho > 0 {
        let rn = if first {
            first = false;
            (nkouho as f64 * 0.5) as usize
        } else {
            (st.rng.rnd() * nkouho as f64) as usize
        };
        let picktmp = pickkouho[rn];
        nkouho -= 1;
        pickkouho[rn] = pickkouho[nkouho];
        if !picks.iter().any(|&p| same_at(st, picktmp, p)) {
            picks.push(picktmp);
        }
    }
    picks.sort_unstable();
    let npick = picks.len();

    // `pickmtx` (with-diagonal half matrix, `:1620-1713`). Row 0 reuses the
    // reference distances of the sorted entries.
    let mut pickmtx: Vec<Vec<f64>> = (0..npick).map(|i| vec![0.0f64; npick - i]).collect();
    for k in 0..npick {
        pickmtx[0][k] = scores[picks[k]].score;
    }
    for j in 1..npick {
        let r = scores[picks[j]].numinseq;
        for i in (j + 1)..npick {
            pickmtx[j][i - j] = st.scorer.dist(r, scores[picks[i]].numinseq);
        }
    }
    let pick_of: Vec<Option<usize>> = {
        let mut m = vec![None; nin];
        for (k, &p) in picks.iter().enumerate() {
            m[p] = Some(k);
        }
        m
    };

    // Redundancy filter (`HUKINTOTREE`, `:1907-1982`).
    let mut tsukau = vec![true; npick];
    if npick > 2 {
        let kijun = maxdist * st.params.tokyoripara;
        for i in 0..npick - 1 {
            if !tsukau[i] { continue; }
            for j in (i + 1)..npick {
                if !tsukau[j] { continue; }
                if pickmtx[i][j - i] < kijun {
                    tsukau[j] = false;
                }
            }
        }
    }
    let yuko_picks: Vec<usize> = (0..npick).filter(|&i| tsukau[i]).collect();
    let nyuko = yuko_picks.len();
    let yukos: Vec<usize> = yuko_picks.iter().map(|&i| picks[i]).collect();

    let mut outs: Vec<Vec<SplitEntry>> = vec![Vec::new(); nyuko];
    let mut treeorder: Vec<usize> = vec![0];
    let mut topol: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    let mut yuko_lengths: Vec<(f64, f64)> = Vec::new();
    let mut child_uniform_leaf = uniform_leaf;

    if nyuko == 1 {
        // Every pick but the reference was a duplicate of it
        // (`:2105-2123`): one group, emitted as a leaf.
        child_uniform_leaf = true;
        outs[0] = scores.clone();
    } else {
        // `dfromc[yuko][j]` (`:1715-1729`, `:2134-2236`). Columns that are
        // picks come from `pickmtx`; the rest are computed.
        let mut dfromc: Vec<Vec<f64>> = vec![vec![0.0f64; nin]; nyuko];
        for j in 0..nin {
            dfromc[0][j] = scores[j].score;
        }
        for (yi, &p) in yuko_picks.iter().enumerate().skip(1) {
            let r = scores[picks[p]].numinseq;
            for j in 0..nin {
                dfromc[yi][j] = match pick_of[j] {
                    Some(k) => {
                        let (a, b) = if p <= k { (p, k) } else { (k, p) };
                        pickmtx[a][b - a]
                    }
                    None => st.scorer.dist(r, scores[j].numinseq),
                };
            }
        }

        // Assignment (`:2242-2303`): yukos belong to themselves; others to
        // an identical yuko if any, else the strictly closest one.
        let mut s_y_map: Vec<Option<usize>> = vec![None; nin];
        for (yi, &y) in yukos.iter().enumerate() {
            s_y_map[y] = Some(yi);
        }
        for j in 0..nin {
            let belong = match s_y_map[j] {
                Some(yi) => yi,
                None => {
                    let mut b = 0usize;
                    let mut minscore = dfromc[0][j];
                    for yi in 0..nyuko {
                        if st.scorer.same(scores[j].numinseq, scores[yukos[yi]].numinseq) {
                            b = yi;
                            break;
                        }
                        if dfromc[yi][j] < minscore {
                            minscore = dfromc[yi][j];
                            b = yi;
                        }
                    }
                    b
                }
            };
            outs[belong].push(scores[j].clone());
        }

        // UPGMA on the yuko distances (`:2311-2326`).
        if nyuko > 2 {
            let mut dm = DistanceMatrix::new(nyuko);
            for (ii, &i) in yuko_picks.iter().enumerate() {
                for (jj, &j) in yuko_picks.iter().enumerate().skip(ii + 1) {
                    dm.set(ii, jj, pickmtx[i][j - i]);
                }
            }
            let yt = musclesupg(&dm, ClusterMethod::Mix { sueff: 0.1 });
            topol = c_topol(&yt.steps);
            yuko_lengths = yt.steps.iter().map(|s| (s.left_length, s.right_length)).collect();
        } else {
            topol = vec![(vec![0], vec![1])];
            yuko_lengths = vec![(0.0, 0.0)];
        }
        let root = &topol[nyuko - 2];
        treeorder = root.0.iter().chain(root.1.iter()).copied().collect();
    }

    // Recurse into each group in `treeorder` (`:2362-2378`).
    let child_qinoya = Some(scores[0].numinseq);
    let mut children: Vec<String> = vec![String::new(); nyuko];
    for &ii in treeorder.iter().take(nyuko) {
        let group = std::mem::take(&mut outs[ii]);
        let members: Vec<usize> = group.iter().map(|s| s.numinseq).collect();
        children[ii] = splitseq_mq(st, group, child_uniform_leaf, child_qinoya);
        outs[ii] = members
            .into_iter()
            .map(|m| SplitEntry { numinseq: m, selfscore: 0, orilen: 0, score: 0.0 })
            .collect();
    }

    if nyuko == 1 {
        return std::mem::take(&mut children[0]);
    }

    // Merge the yuko groups along the UPGMA tree (`:2444-2554`).
    let mut parttree: Vec<Option<String>> = children.into_iter().map(Some).collect();
    let mut tree = String::new();
    for (l, (t0, t1)) in topol.iter().enumerate() {
        let collect = |side: &[usize]| -> Vec<usize> {
            let mut m: Vec<usize> = side
                .iter()
                .flat_map(|&y| outs[y].iter().map(|s| s.numinseq))
                .collect();
            m.sort_unstable();
            m
        };
        let mem1 = collect(t0);
        let mem2 = collect(t1);
        let (len1, len2) = yuko_lengths[l];
        let step = if mem1[0] < mem2[0] {
            JoinStep { left: mem1, right: mem2, left_length: len1, right_length: len2 }
        } else {
            JoinStep { left: mem2, right: mem1, left_length: len2, right_length: len1 }
        };
        st.steps.push(step);

        let (mut v1, mut v2) = (t0[0], t1[0]);
        if v1 > v2 {
            std::mem::swap(&mut v1, &mut v2);
        }
        let s1 = parttree[v1].take().unwrap_or_default();
        let s2 = parttree[v2].take().unwrap_or_default();
        tree = format!("({},{})", s1, s2);
        parttree[v1] = Some(tree.clone());
    }
    tree
}

/// Length correction of `splittbfast.c:1409` (and `:1674`, `:2170`).
/// arm64 clang contracts `shorter / longer * lenfacd + lenfacb / ( longer
/// + lenfacc )` into one fused multiply-add, so it goes through
/// [`mafft_types::fp::fmadd`].
fn lenfac_splittbfast(len1: usize, len2: usize, a: f64, b: f64, c: f64, d: f64) -> f64 {
    let (longer, shorter) = if len1 > len2 {
        (len1 as f64, len2 as f64)
    } else {
        (len2 as f64, len1 as f64)
    };
    1.0 / (mafft_types::fp::fmadd(shorter / longer, d, b / (longer + c)) + a)
}

/// `splittbfast` default (k-tuple) distance: `localcommonsextet_p` on the
/// 6-mer point tables with the `lenfac` correction, clamped at `MAX6DIST`.
pub struct KtupleScorer<'a> {
    pub sequences: &'a [Vec<u8>],
    pub points: Vec<Vec<u32>>,
    pub selfscore: Vec<i64>,
    pub tsize: usize,
    pub lenfac_abcd: (f64, f64, f64, f64),
    table_cache: RefCell<Option<(usize, Vec<i32>)>>,
}

impl<'a> KtupleScorer<'a> {
    pub fn new(sequences: &'a [Vec<u8>], kind: crate::parttree_pivot::PtSeqKind) -> Self {
        use crate::parttree_dist::*;
        use crate::parttree_pivot::PtSeqKind;
        let (tsize, abcd) = match kind {
            PtSeqKind::Protein => (46656, (PLENFACA, PLENFACB, PLENFACC, PLENFACD)),
            PtSeqKind::Dna => (4096, (DLENFACA, DLENFACB, DLENFACC, DLENFACD)),
        };
        let points: Vec<Vec<u32>> = sequences
            .iter()
            .map(|s| match kind {
                PtSeqKind::Protein => encode_points_protein(s),
                PtSeqKind::Dna => encode_points_dna(s),
            })
            .collect();
        let selfscore = points.iter().map(|p| p.len() as i64).collect();
        Self {
            sequences,
            points,
            selfscore,
            tsize,
            lenfac_abcd: abcd,
            table_cache: RefCell::new(None),
        }
    }

    pub fn entries(&self) -> Vec<SplitEntry> {
        (0..self.sequences.len())
            .map(|i| SplitEntry {
                numinseq: i,
                selfscore: self.selfscore[i],
                orilen: self.sequences[i].len(),
                score: 0.0,
            })
            .collect()
    }
}

impl SplitScorer for KtupleScorer<'_> {
    fn dist(&self, r: usize, o: usize) -> f64 {
        let mut cache = self.table_cache.borrow_mut();
        if cache.as_ref().map(|(k, _)| *k) != Some(r) {
            *cache = Some((r, composition_table(&self.points[r], self.tsize)));
        }
        let table = &cache.as_ref().unwrap().1;
        let common = common_sextets_p(table, &self.points[o], self.tsize);
        let bunbo = self.selfscore[r].min(self.selfscore[o]) as f64;
        let (a, b, c, d) = self.lenfac_abcd;
        let lf = lenfac_splittbfast(self.sequences[r].len(), self.sequences[o].len(), a, b, c, d);
        let raw = if bunbo > 0.0 { 1.0 - common as f64 / bunbo } else { 1.0 };
        let dist = raw * lf;
        if dist > MAX6DIST { MAX6DIST } else { dist }
    }

    fn same(&self, a: usize, b: usize) -> bool {
        self.sequences[a] == self.sequences[b]
    }
}

/// `-Z` (`fromaln`) distance: `naivepairscore11` on the input alignment
/// (`splittbfast.c:1428`, `:1693`, `:2213`); no clamping.
pub struct FromalnScorer<'a> {
    pub aligned: &'a [Vec<u8>],
    pub degapped: Vec<Vec<u8>>,
    pub selfscore: Vec<i64>,
    pub matrix: &'a [Vec<f64>],
    pub amino_map: &'a [u8; 256],
    pub penalty: f64,
}

impl<'a> FromalnScorer<'a> {
    pub fn new(
        aligned: &'a [Vec<u8>],
        matrix: &'a [Vec<f64>],
        amino_map: &'a [u8; 256],
        penalty: f64,
    ) -> Self {
        let degapped: Vec<Vec<u8>> = aligned
            .iter()
            .map(|s| s.iter().copied().filter(|&c| c != b'-').collect())
            .collect();
        let nalpha = matrix.len();
        // `splittbfast.c:3000-3006`: sum of `amino_dis[c][c]`.
        let selfscore = degapped
            .iter()
            .map(|s| {
                let mut v = 0.0f64;
                for &c in s {
                    let i = amino_map[c as usize] as usize;
                    if i < nalpha {
                        v += matrix[i][i];
                    }
                }
                v as i64
            })
            .collect();
        Self { aligned, degapped, selfscore, matrix, amino_map, penalty }
    }

    pub fn entries(&self) -> Vec<SplitEntry> {
        (0..self.aligned.len())
            .map(|i| SplitEntry {
                numinseq: i,
                selfscore: self.selfscore[i],
                orilen: self.degapped[i].len(),
                score: 0.0,
            })
            .collect()
    }
}

/// `mltaln9.c::naivepairscore11`: common gaps removed, one `penal` per gap
/// run, `amino_dis` elsewhere.
pub fn naivepairscore11(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penal: f64,
) -> f64 {
    let (s1, s2): (Vec<u8>, Vec<u8>) = seq1
        .iter()
        .zip(seq2.iter())
        .filter(|&(&a, &b)| !(a == b'-' && b == b'-'))
        .map(|(&a, &b)| (a, b))
        .unzip();
    let n = s1.len();
    let nalpha = matrix.len();
    let mut vali = 0.0f64;
    let mut k = 0usize;
    while k < n {
        if s1[k] == b'-' {
            vali += penal;
            while k < n && s1[k] == b'-' {
                k += 1;
            }
            continue;
        }
        if s2[k] == b'-' {
            vali += penal;
            while k < n && s2[k] == b'-' {
                k += 1;
            }
            continue;
        }
        let i = amino_map[s1[k] as usize] as usize;
        let j = amino_map[s2[k] as usize] as usize;
        if i < nalpha && j < nalpha {
            vali += matrix[i][j];
        }
        k += 1;
    }
    vali
}

impl SplitScorer for FromalnScorer<'_> {
    fn dist(&self, r: usize, o: usize) -> f64 {
        let pair = naivepairscore11(
            &self.aligned[o], &self.aligned[r], self.matrix, self.amino_map, self.penalty,
        );
        let bunbo = self.selfscore[r].min(self.selfscore[o]) as f64;
        if bunbo > 0.0 { 1.0 - pair / bunbo } else { 1.0 }
    }

    fn same(&self, a: usize, b: usize) -> bool {
        self.degapped[a] == self.degapped[b]
    }
}

/// Closure-backed scorer (used for the `-U` DP distance of `--dpparttree`).
pub struct FnScorer<D, E>
where
    D: Fn(usize, usize) -> f64,
    E: Fn(usize, usize) -> bool,
{
    pub dist: D,
    pub same: E,
}

impl<D, E> SplitScorer for FnScorer<D, E>
where
    D: Fn(usize, usize) -> f64,
    E: Fn(usize, usize) -> bool,
{
    fn dist(&self, r: usize, o: usize) -> f64 {
        (self.dist)(r, o)
    }
    fn same(&self, a: usize, b: usize) -> bool {
        (self.same)(a, b)
    }
}

/// `--parttree` CALL 1: k-tuple `splittbfast` on the raw sequences.
pub fn run_ktuple(sequences: &[Vec<u8>], kind: crate::parttree_pivot::PtSeqKind) -> SplitRun {
    let scorer = KtupleScorer::new(sequences, kind);
    let entries = scorer.entries();
    run_splitseq(entries, &scorer, SplitParams::for_njob(sequences.len(), false))
}

/// `--parttree` CALL 2+: `splittbfast -Z` on the previous call's alignment.
pub fn run_fromaln(
    aligned: &[Vec<u8>],
    matrix: &[Vec<f64>],
    amino_map: &[u8; 256],
    penalty: f64,
) -> SplitRun {
    let scorer = FromalnScorer::new(aligned, matrix, amino_map, penalty);
    let entries = scorer.entries();
    run_splitseq(entries, &scorer, SplitParams::for_njob(aligned.len(), true))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bsd_rand_matches_macos_libc() {
        let mut r = BsdRand::default();
        let got: Vec<u32> = (0..5).map(|_| r.rand()).collect();
        assert_eq!(got, vec![16807, 282475249, 1622650073, 984943658, 1144108930]);
    }
}
