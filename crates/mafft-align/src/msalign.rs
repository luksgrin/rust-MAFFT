//! Hirschberg-style linear-space profile DP — port of C's `MSalignmm`
//! (`mafft-upstream/core/MSalignmm.c`).
//!
//! C MAFFT routes `--memsave` (`alg = 'M'`, `tbfast.c:1159-1161`) and
//! the auto-switch path at `len > 30000` (`tbfast.c:1096-1104`) through
//! `MSalignmm`. The DP recurrence is identical to `A__align` /
//! [`profile_align_imp_with_boundary`]; the only difference is memory
//! layout — `A__align` allocates an O(N×M) trace matrix, `MSalignmm`
//! uses Hirschberg's divide-and-conquer to work in O(N+M).
//!
//! Implementation strategy:
//!   * Base case (`lgth1 < DPTANNI=100` OR `lgth2 < DPTANNI`): delegate
//!     to [`profile_align_imp_with_boundary`]. The full DP is fine for
//!     small sub-problems.
//!   * Recursive case: forward DP from row `ist` to row `imid =
//!     ist + lgth1/2`, capturing `midw[j]`, `midm[j]`, `midn[j]` for
//!     each column j at the midpoint row. Backward DP from row `ien`
//!     to `imid`, adding to those arrays. Find `jmid` = argmax over
//!     midw/midm/midn. Recurse on top-left `(ist, ist+jumpi, jst,
//!     jst+jumpj)` and bottom-right `(ist+imid, ien, jst+jmid, jen)`
//!     with appropriate inter-half gap padding.
//!
//! ## Correctness scope
//!
//! Verified:
//!   * Base case (< DPTANNI on either axis) is byte-identical to
//!     [`profile_align_imp_with_boundary`], because the recursion
//!     bottoms out directly into it.
//!   * Recursive case returns the SAME optimal score as the full DP
//!     for every test input — both for identical-input traces (no
//!     ties) and for gap-required traces with unique optima.
//!   * For inputs with a unique optimal alignment, the trace is
//!     byte-identical to the full DP.
//!
//! An m-split (`jumpi = jumpbackj[jmid]` above the midpoint) continues
//! the backward DP down to that row before reading `jumpforwi` /
//! `jumpforwj`, and a subproblem with a non-positive length on one
//! side emits the other side as terminal gaps (`MSalignmm_rec` when
//! `lgth2 <= 0`, and `Atracking` when `lgth1 <= 0`). Stopping at the
//! midpoint, or skipping that call, returns an alignment narrower
//! than the input profiles.

use mafft_types::fp::fmadd;
use crate::dp::{Alignment, AlignOp, GapModel};
use crate::profile::{Profile, profile_align_imp_with_boundary, BoundaryFreqs, build_sparse_cpmx};

/// C `MSalignmm.c:11` — sub-problem size threshold below which the
/// recursion bottoms out into the direct DP.
const DPTANNI: usize = 100;

/// Top-level entry point. Mirrors C `MSalignmm` (line 2094): set up
/// `headgapfreq1/2` defaults (= 1.0, since the `sgap1`/`egap1` path
/// isn't relevant for `--memsave`), then drive `msalignmm_rec`.
pub fn msalignmm(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<f64>],
    gap: &GapModel,
    head_gap: bool,
    tail_gap: bool,
) -> Alignment {
    msalignmm_with_boundary(prof1, prof2, matrix, gap, head_gap, tail_gap, BoundaryFreqs::default())
}

/// [`msalignmm`] on a cut segment. `MSalignmm` with `sgap`/`egap`
/// (`MSalignmm.c:2331-2341`, called from `Falign_udpari_long`,
/// `Falign.c:3354`) takes `headgapfreq` from the column before the cut
/// and `gapfreq[lgth]` from the column after it, as nongap fractions.
pub fn msalignmm_with_boundary(
    prof1: &Profile,
    prof2: &Profile,
    matrix: &[Vec<f64>],
    gap: &GapModel,
    head_gap: bool,
    tail_gap: bool,
    boundary: BoundaryFreqs,
) -> Alignment {
    let n = prof1.length;
    let m = prof2.length;
    if n == 0 || m == 0 {
        return Alignment::against_gaps(n, m);
    }
    let mut ctx = MsCtx::new(prof1, prof2, matrix, gap, boundary.tail1, boundary.tail2);
    let mut ops = Vec::with_capacity(n + m);
    let score = msalignmm_rec(
        &mut ctx,
        0, n - 1, 0, m - 1,
        head_gap, tail_gap,
        boundary.head1, boundary.head2,
        &mut ops,
    );
    debug_assert_eq!(
        ops.iter().filter(|op| matches!(op, AlignOp::Match | AlignOp::Delete)).count(),
        n,
        "msalignmm dropped columns of profile 1",
    );
    debug_assert_eq!(
        ops.iter().filter(|op| matches!(op, AlignOp::Match | AlignOp::Insert)).count(),
        m,
        "msalignmm dropped columns of profile 2",
    );
    Alignment {
        seq1: Vec::new(),
        seq2: Vec::new(),
        score,
        operations: ops,
    }
}

/// Per-call state shared by every node of the recursion.
///
/// C `MSalignmm` computes the gap-cost vectors and the nongap fractions
/// once per call (`MSalignmm.c:2298-2311`) and passes them down as
/// `gapinfo`; every DP row then goes through
/// `match_calc_alphabet_seq` (`MSalignmm.c:182-236`), which visits only
/// the non-zero alphabet entries of each profile column
/// (`cpmxpd` / `cpmxpdn`, built once per direction with
/// `initialize = 1`). `sp1` / `sp2` are those sparse views over the
/// whole of `prof1` / `prof2`; a node reads its own window through
/// the offsets table.
struct MsCtx<'a> {
    prof1: &'a Profile,
    prof2: &'a Profile,
    matrix: &'a [Vec<f64>],
    gap: &'a GapModel,
    nalpha: usize,
    /// `gapinfo[0..4]` in absolute coordinates:
    /// `0.5 * (1 - count) * penalty * nongap`, plus one zero of calloc
    /// padding at `[length]` (C reads `ogcp[lgth]` on the last column).
    ogcp1: Vec<f64>,
    fgcp1: Vec<f64>,
    ogcp2: Vec<f64>,
    fgcp2: Vec<f64>,
    /// `gapfreq{1,2}f`: the nongap fraction of each column, and at
    /// `[length]` the segment's tail value (`outgapcount` of `egap`,
    /// 1.0 for a whole profile).
    gf1: Vec<f64>,
    gf2: Vec<f64>,
    /// Sparse column views (C `cpmxpd` / `cpmxpdn`): column `j`'s
    /// non-zero `(alphabet_index, freq)` pairs, in alphabet order, are
    /// `ent[off[j]..off[j + 1]]`.
    sp1_off: Vec<usize>,
    sp1_ent: Vec<(usize, f64)>,
    sp2_off: Vec<usize>,
    sp2_ent: Vec<(usize, f64)>,
    /// `scarr` of `match_calc_alphabet_seq`, reused across rows.
    scarr: Vec<f64>,
}

impl<'a> MsCtx<'a> {
    fn new(
        prof1: &'a Profile,
        prof2: &'a Profile,
        matrix: &'a [Vec<f64>],
        gap: &'a GapModel,
        tail1: f64,
        tail2: f64,
    ) -> Self {
        let nalpha = prof1.nalphabets.min(prof2.nalphabets).min(matrix.len());
        // C `MSalignmm.c:2298-2306`:
        //   ogcp1[i] = 0.5 * ( 1.0 - ogcp1opt[i] ) * fpenalty * gapfreq1pt[i]
        // `Profile.ogcp` / `.fgcp` hold the raw opening / closing counts.
        let effective = |count: &[f64], p: &Profile| -> Vec<f64> {
            let mut v: Vec<f64> = (0..p.length)
                .map(|i| 0.5 * (1.0 - count[i]) * gap.open * p.nongap_freq[i])
                .collect();
            v.push(0.0);
            v
        };
        let ogcp1 = effective(&prof1.ogcp, prof1);
        let fgcp1 = effective(&prof1.fgcp, prof1);
        let ogcp2 = effective(&prof2.ogcp, prof2);
        let fgcp2 = effective(&prof2.fgcp, prof2);
        let mut gf1 = prof1.nongap_freq[..prof1.length].to_vec();
        gf1.push(tail1);
        let mut gf2 = prof2.nongap_freq[..prof2.length].to_vec();
        gf2.push(tail2);
        let (mut sp1_off, mut sp1_ent) = (Vec::new(), Vec::new());
        let (mut sp2_off, mut sp2_ent) = (Vec::new(), Vec::new());
        build_sparse_cpmx(&prof1.freqs, prof1.length, nalpha, &mut sp1_off, &mut sp1_ent);
        build_sparse_cpmx(&prof2.freqs, prof2.length, nalpha, &mut sp2_off, &mut sp2_ent);
        Self {
            prof1, prof2, matrix, gap, nalpha,
            ogcp1, fgcp1, ogcp2, fgcp2,
            gf1, gf2,
            sp1_off, sp1_ent, sp2_off, sp2_ent,
            scarr: vec![0.0; nalpha],
        }
    }

    /// `match_calc_alphabet_seq( mtx, out, cpmx1pt, cpmx2pt, row, jst,
    /// lgth2, .. )`: `out[j]` is the profile match score of `prof1`
    /// column `row` against `prof2` column `jst + j`.
    fn match_row(&mut self, row: usize, jst: usize, lgth2: usize, out: &mut [f64]) {
        let nalpha = self.nalpha;
        match_into(
            &mut self.scarr[..nalpha], self.matrix, &self.prof1.freqs[row][..nalpha],
            &self.sp2_off, &self.sp2_ent, jst, lgth2, out,
        );
    }

    /// The transposed call (`cpmx2pt, cpmx1pt, col, ist, lgth1`):
    /// `out[i]` is `prof2` column `col` against `prof1` column `ist + i`.
    fn match_col(&mut self, col: usize, ist: usize, lgth1: usize, out: &mut [f64]) {
        let nalpha = self.nalpha;
        match_into(
            &mut self.scarr[..nalpha], self.matrix, &self.prof2.freqs[col][..nalpha],
            &self.sp1_off, &self.sp1_ent, ist, lgth1, out,
        );
    }
}

/// Body of C `match_calc_alphabet_seq` (`MSalignmm.c:216-235`).
///
/// `scarr[l] = sum_k matrix[k][l] * freqs1[k]`, accumulated over `k` in
/// ascending order for every `l` (the nest is k-outer / l-inner so the
/// matrix rows are walked contiguously; each `scarr[l]` still sees the
/// same ordered multiply-add sequence). Then `out[j]` is the sum, in
/// alphabet order, of `scarr[k] * v` over column `start2 + j`'s non-zero
/// entries only. A zero entry would contribute `fma(scarr, 0, acc)`,
/// which is `acc` exactly, so skipping it leaves the bit pattern
/// unchanged relative to the dense sum.
#[allow(clippy::too_many_arguments)]
#[inline]
fn match_into(
    scarr: &mut [f64],
    matrix: &[Vec<f64>],
    freqs1: &[f64],
    off: &[usize],
    ent: &[(usize, f64)],
    start2: usize,
    lgth2: usize,
    out: &mut [f64],
) {
    let nalpha = scarr.len();
    scarr.fill(0.0);
    for (mrow, &f1) in matrix[..nalpha].iter().zip(freqs1) {
        for (s, &mv) in scarr.iter_mut().zip(&mrow[..nalpha]) {
            *s = fmadd(mv, f1, *s);
        }
    }
    let offs = &off[start2..=start2 + lgth2];
    let mut ents = ent[offs[0]..offs[lgth2]].iter();
    for (o, w) in out[..lgth2].iter_mut().zip(offs.windows(2)) {
        let mut acc = 0.0f64;
        for &(k, v) in ents.by_ref().take(w[1] - w[0]) {
            acc = fmadd(scarr[k], v, acc);
        }
        *o = acc;
    }
}

/// Recursive Hirschberg DP. `ist`/`ien` and `jst`/`jen` are inclusive
/// absolute positions into the parent `prof1` / `prof2`. The function
/// appends the trace operations for the sub-region into `out_ops` and
/// returns the score for that sub-region.
///
/// Mirrors C `MSalignmm_rec` (line 991). The `headgapfreq{1,2}_g`
/// parameters carry the "nongap fraction at position ist-1 / jst-1"
/// when those positions don't exist in the parent profile (i.e. for
/// the top-level call); for nested calls the actual previous column's
/// nongap_freq is used (matching C line 1072-1075).
#[allow(clippy::too_many_arguments)]
fn msalignmm_rec(
    ctx: &mut MsCtx<'_>,
    ist: usize, ien: usize,
    jst: usize, jen: usize,
    head_gap: bool, tail_gap: bool,
    headgapfreq1_g: f64, headgapfreq2_g: f64,
    out_ops: &mut Vec<AlignOp>,
) -> f64 {
    // C `MSalignmm_rec:1099`: `lgth2 <= 0` copies the remaining seq1
    // columns and fills seq2 with gaps. A non-positive `lgth1` falls
    // through to `MSalignmm_tanni` / `Atracking`, which emits the
    // remaining seq2 columns as gaps in seq1 (`ijp[0][j] = -(j+1)`).
    // Skipping either call drops those columns, so the alignment is
    // narrower than its inputs.
    if ist > ien || jst > jen {
        if ist <= ien {
            for _ in ist..=ien {
                out_ops.push(AlignOp::Delete);
            }
        } else if jst <= jen {
            for _ in jst..=jen {
                out_ops.push(AlignOp::Insert);
            }
        }
        return 0.0;
    }
    let lgth1 = ien - ist + 1;
    let lgth2 = jen - jst + 1;

    // Base case: small sub-problem → full DP.
    // C `MSalignmm_rec:1144` — `if( lgth1 < DPTANNI || lgth2 < DPTANNI )`.
    if lgth1 < DPTANNI || lgth2 < DPTANNI {
        return base_case(
            ctx,
            ist, ien, jst, jen,
            head_gap, tail_gap,
            headgapfreq1_g, headgapfreq2_g,
            out_ops,
        );
    }

    // C `MSalignmm_rec:1261` — `imid = lgth1 * 0.5;` (integer truncation).
    let imid = lgth1 / 2;

    // Forward DP rows 1..=imid, accumulating midw/midm/midn at row imid.
    let fwd = forward_dp(ctx, ist, jst, jen, imid, headgapfreq1_g, headgapfreq2_g);

    // Backward DP rows lgth1-2 down to imid-1, adding to midw/midm/midn.
    let bwd = backward_dp(ctx, ist, ien, jst, jen, imid, fwd);

    let (jumpi, jumpj, jmid) = bwd.split_point;
    let imid = bwd.imid_override.unwrap_or(imid);

    // Recurse on top-left sub-region `(ist .. ist+jumpi, jst .. jst+jumpj)`.
    // C `MSalignmm_rec:1886` — both endpoints inclusive. For a diagonal
    // join (`midw` winner), `jumpi = imid-1` and `jumpj = jmid-1`, so
    // top covers `[ist, ist+imid-1]` (imid rows). Bottom-right starts
    // at `ist+imid` so the union is exact.
    // C passes the caller's `headgp` and `tailgp` to both halves
    // (`MSalignmm.c:1886`, `:1975`); the base case decides from the
    // window position whether a boundary is terminal.
    let mut value = 0.0;
    value += msalignmm_rec(
        ctx,
        ist, ist + jumpi, jst, jst + jumpj,
        head_gap, tail_gap,
        headgapfreq1_g, headgapfreq2_g,
        out_ops,
    );

    // Inter-half gap padding for prof1: `l = jmid - jumpj - 1`
    // (C line 1907). This is a run of gap-in-prof1 / residue-in-prof2.
    let l_horiz: i64 = jmid as i64 - jumpj as i64 - 1;
    if l_horiz > 0 {
        for _ in 0..l_horiz {
            out_ops.push(AlignOp::Insert);
        }
        // Gap-cost contribution (C line 1925):
        //   value += ogcp2[jumpj+1] + fgcp2[jmid-1]
        // In C, `ogcp2 = gapinfo[2] + jst` is a SLICE starting at jst,
        // so `ogcp2[jumpj+1]` is `ogcp2_full[jst + jumpj + 1]` in
        // absolute coordinates. Same for fgcp2.
        let idx_o = jst + jumpj + 1;
        let idx_f = jst + jmid - 1;
        if idx_o < ctx.ogcp2.len() && idx_f < ctx.fgcp2.len() {
            value += ctx.ogcp2[idx_o] + ctx.fgcp2[idx_f];
        }
    }
    // Inter-half gap padding for prof2: `l = imid - jumpi - 1`
    // (C line 1929).
    let l_vert: i64 = imid as i64 - jumpi as i64 - 1;
    if l_vert > 0 {
        for _ in 0..l_vert {
            out_ops.push(AlignOp::Delete);
        }
        // C line 1953:  value += ogcp1[jumpi+1] + fgcp1[imid-1]
        // — slice-relative, same +1 / -1 convention as prof2.
        let idx_o = ist + jumpi + 1;
        let idx_f = ist + imid - 1;
        if idx_o < ctx.ogcp1.len() && idx_f < ctx.fgcp1.len() {
            value += ctx.ogcp1[idx_o] + ctx.fgcp1[idx_f];
        }
    }

    // Recurse on bottom-right sub-region `(ist+imid .. ien, jst+jmid .. jen)`.
    // For the bottom recursion the head boundary lives INSIDE the
    // parent: C `MSalignmm_tanni` derives `headgapfreq1 = gapfreq1f[-1]`
    // from the original profile when `ist > 0` (line 709-712). We
    // mirror that here by passing the previous-column nongap_freq.
    // Always recurse. A start past the parent's end is the `lgth <= 0`
    // case above (C still calls `_rec`; the empty side becomes gaps).
    // Gating on both bounds being in range dropped the other profile's
    // remaining columns whenever `jmid == lgth2`.
    let bottom_ist = ist + imid;
    let bottom_jst = jst + jmid;
    let head1_for_bot = if bottom_ist > 0 && bottom_ist - 1 < ctx.prof1.length {
        ctx.gf1[bottom_ist - 1]
    } else { headgapfreq1_g };
    let head2_for_bot = if bottom_jst > 0 && bottom_jst - 1 < ctx.prof2.length {
        ctx.gf2[bottom_jst - 1]
    } else { headgapfreq2_g };
    value += msalignmm_rec(
        ctx,
        bottom_ist, ien, bottom_jst, jen,
        head_gap, tail_gap,
        head1_for_bot, head2_for_bot,
        out_ops,
    );

    value
}

/// Base case: sub-region small enough that the full O(lgth1 × lgth2)
/// trace matrix fits comfortably. Use the existing
/// [`profile_align_imp_with_boundary`] (= our port of C `A__align` /
/// `MSalignmm_tanni` — they share the same recurrence).
#[allow(clippy::too_many_arguments)]
fn base_case(
    ctx: &mut MsCtx<'_>,
    ist: usize, ien: usize,
    jst: usize, jen: usize,
    head_gap: bool, tail_gap: bool,
    headgapfreq1_g: f64, headgapfreq2_g: f64,
    out_ops: &mut Vec<AlignOp>,
) -> f64 {
    let sub1 = ctx.prof1.sub_profile(ist, ien + 1);
    let sub2 = ctx.prof2.sub_profile(jst, jen + 1);

    // Derive boundary nongap fractions from the parent profile.
    // C `MSalignmm_tanni:709-712`:
    //   if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
    //   else headgapfreq1 = headgapfreq1_g;
    // (and similarly for tail). `gapfreq1f` is `1 - gap_freq` per
    // position (the "nongap fraction"); `gapfreq1f[-1]` peeks at
    // position `ist - 1` in the parent profile.
    let head1 = if ist > 0 { ctx.gf1[ist - 1] } else { headgapfreq1_g };
    let head2 = if jst > 0 { ctx.gf2[jst - 1] } else { headgapfreq2_g };
    let tail1 = ctx.gf1[ien + 1];
    let tail2 = ctx.gf2[jen + 1];

    // Head penalties follow `MSalignmm_tanni:781-788`: an internal
    // start (`ist != 0` or `jst != 0`) is charged even when `headgp`
    // is 0. Tail freedom follows `Atracking` (`MSalignmm.c:512-514`):
    // with `tailgp == 1` the tail is always charged; with `tailgp == 0`
    // the free-end scan runs when EITHER window reaches the global end
    // (`ien == fulllen1-1 || jen == fulllen2-1`). Requiring both
    // penalized the tail of a fragment that ends inside a longer
    // profile, which slid its last residue to a later exact match.
    let effective_head = head_gap || ist != 0 || jst != 0;
    let at_global_end = ien + 1 == ctx.prof1.length || jen + 1 == ctx.prof2.length;
    let effective_tail = tail_gap || !at_global_end;

    let aln = profile_align_imp_with_boundary(
        &sub1, &sub2, ctx.matrix, ctx.gap, effective_head, effective_tail, None, false,
        BoundaryFreqs { head1, head2, tail1, tail2 },
    );
    out_ops.extend(aln.operations);
    aln.score
}

/// Result of the forward DP up to and including row `imid`.
struct ForwardState {
    /// `midw[j]` after the forward pass at row `imid`:
    /// `currentw[j] + match_calc` accumulated up to row `imid`.
    midw: Vec<f64>,
    /// `midm[j]` after the forward pass: the row-running `m[j]` value
    /// (gap-extend tracker in prof2 direction) at row `imid`.
    midm: Vec<f64>,
    /// `midn[j]` after the forward pass: the col-running `mi` value
    /// (gap-extend tracker in prof1 direction) at row `imid`.
    midn: Vec<f64>,
    /// Trace metadata for the forward sweep at row `imid`:
    /// `jumpbackj[j] = *mpjpt` at the column update (i.e. the prior i
    /// where mj[j] was last raised), `jumpbacki[j] = mpi` (the prior j
    /// where mi was last raised).
    jumpbackj: Vec<i64>,
    jumpbacki: Vec<i64>,
    /// Row buffers handed on to the backward pass for reuse.
    currentw: Vec<f64>,
    previousw: Vec<f64>,
    m: Vec<f64>,
    mp: Vec<i64>,
}

/// Result of the backward DP plus the split-point selection.
struct BackwardState {
    /// `(jumpi, jumpj, jmid)` — the optimal split coordinates.
    /// `jumpi`/`jumpj` are the i/j right BEFORE the cut (top-left
    /// half ends at `ist+jumpi, jst+jumpj` inclusive); `imid`/`jmid`
    /// are the i/j right AFTER the cut (bottom-right half starts at
    /// `ist+imid, jst+jmid`).
    split_point: (usize, usize, usize),
    /// If `Some`, overrides the parent's `imid = lgth1/2` with
    /// `jumpforwi[jumpj]` (C line 1767). For typical splits this
    /// equals `imid` (no change); for splits where the mj-branch
    /// fed `ijpi = mp[dj]` it points further down, shifting the
    /// bottom-half recursion start row.
    imid_override: Option<usize>,
}

/// Forward DP from row `ist` to row `ist + imid` (inclusive of the
/// midpoint sweep). Mirrors C `MSalignmm_rec:1220-1413` — the same
/// recurrence as `MSalignmm_tanni` / [`profile_align_imp_with_boundary`]
/// but without the O(N×M) trace matrix.
///
/// Every row index touched here is at most `ien < prof1.length` and
/// every column index at most `jen + 1 <= prof2.length`, so the
/// absolute slices below are in range; the `[length]` entries of the
/// gap vectors are C's calloc padding and the segment tail value.
fn forward_dp(
    ctx: &mut MsCtx<'_>,
    ist: usize,
    jst: usize, jen: usize,
    imid: usize,
    headgapfreq1_g: f64, headgapfreq2_g: f64,
) -> ForwardState {
    let lgth2 = jen - jst + 1;
    debug_assert!(ist + imid < ctx.prof1.length);
    debug_assert!(jen < ctx.prof2.length);

    // Use the parent profile's nongap_freq for `gapfreq{1,2}f` and the
    // headgapfreq{1,2}_g for the absolute boundary (matching C
    // MSalignmm_rec:1072-1075).
    let headgapfreq1 = if ist > 0 { ctx.gf1[ist - 1] } else { headgapfreq1_g };
    let headgapfreq2 = if jst > 0 { ctx.gf2[jst - 1] } else { headgapfreq2_g };

    // C MSalignmm has USE_PENALTY_EX = 0; there is no per-cell
    // extension term anywhere in this function.

    // `initverticalw`: `match_calc(prof2[jst], prof1[ist..])` then add
    // `ogcp1[ist] * headgapfreq2 + fgcp1[ist+i-1] * nongap_freq2[jst]`.
    // We collapse "size lgth1+1" to "absolute indices ist..=ien" — the
    // first element corresponds to row 0 of the sub-region.
    let mut initverticalw = vec![0.0f64; imid + 2];
    ctx.match_col(jst, ist, imid + 1, &mut initverticalw[..=imid]);
    {
        // `MSalignmm_rec` adds this penalty on every call
        // (`MSalignmm.c:1229-1237`), including the free-terminal-gap
        // case `headgp == 0 && ist == 0`, so the forward pass takes no
        // `headgp`. Only `MSalignmm_tanni` gates on `headgp || ist`
        // (`:781`).
        let gf2_0 = ctx.gf2[jst];
        let ogcp1_0 = ctx.ogcp1[ist];
        let fgcp1 = &ctx.fgcp1[ist..ist + imid];
        for (di, iv) in initverticalw[1..=imid].iter_mut().enumerate() {
            *iv = fmadd(fgcp1[di], gf2_0, fmadd(ogcp1_0, headgapfreq2, *iv));
        }
    }

    // `currentw`: `match_calc(prof1[ist], prof2[jst..])` then add
    // head-gap penalty for prof2 axis.
    let mut currentw = vec![0.0f64; lgth2 + 1];
    let mut previousw = vec![0.0f64; lgth2 + 1];
    ctx.match_row(ist, jst, lgth2, &mut currentw[..lgth2]);
    {
        // Same unconditional init (`MSalignmm.c:1234-1237`).
        let gf1_0 = ctx.gf1[ist];
        let ogcp2_0 = ctx.ogcp2[jst];
        let fgcp2 = &ctx.fgcp2[jst..=jen];
        for (dj, cw) in currentw[1..=lgth2].iter_mut().enumerate() {
            *cw = fmadd(fgcp2[dj], gf1_0, fmadd(ogcp2_0, headgapfreq1, *cw));
        }
    }

    // `m[j] = currentw[j-1] + ogcp1[ist+1] * nongap_freq2[jst+j-1]`
    // C line 1252-1257.
    let mut m = vec![0.0f64; lgth2 + 1];
    let mut mp = vec![0i64; lgth2 + 1];
    {
        let ogcp1_1 = ctx.ogcp1[ist + 1];
        let gf2 = &ctx.gf2[jst..=jen];
        for dj in 1..=lgth2 {
            m[dj] = fmadd(ogcp1_1, gf2[dj - 1], currentw[dj - 1]);
        }
    }

    // Storage for the midpoint row state.
    let mut midw = vec![0.0f64; lgth2 + 1];
    let mut midm = vec![0.0f64; lgth2 + 1];
    let mut midn = vec![0.0f64; lgth2 + 1];
    let mut jumpbackj = vec![0i64; lgth2 + 1];
    let mut jumpbacki = vec![0i64; lgth2 + 1];

    // Forward DP rows 1..=imid (relative to ist).
    // Mirrors C `MSalignmm_rec` inner loop (lines 1268-1413).
    for di in 1..=imid {
        std::mem::swap(&mut previousw, &mut currentw);
        previousw[0] = initverticalw[di - 1];

        let row_i = ist + di;
        let row_im1 = ist + di - 1;
        ctx.match_row(row_i, jst, lgth2, &mut currentw[..lgth2]);
        currentw[0] = initverticalw[di];
        // C line 1301: `m[0] = ogcp1[i]` — for absolute index ist+di.
        m[0] = ctx.ogcp1[row_i];

        // Row-invariant gap terms.
        let gf1_i = ctx.gf1[row_i];
        let gf1_im1 = ctx.gf1[row_im1];
        let ogcp1_i = ctx.ogcp1[row_i];
        let fgcp1_im1 = ctx.fgcp1[row_im1];
        // Column windows, indexed by the local `dj`.
        let ogcp2 = &ctx.ogcp2[jst..=jen + 1];
        let fgcp2 = &ctx.fgcp2[jst..=jen];
        let gf2 = &ctx.gf2[jst..=jen + 1];
        let prev = &previousw[..=lgth2];
        let cur = &mut currentw[..=lgth2];
        let m = &mut m[..=lgth2];
        let mp = &mut mp[..=lgth2];

        // `mi = previousw[0] + ogcp2[jst+1] * nongap_freq1[ist+di-1]`
        let mut mi = fmadd(ogcp2[1], gf1_im1, prev[0]);
        let mut mpi: i64 = 0;

        let capture = di == imid;
        // Capture m[0] at midpoint
        if capture { midm[0] = m[0]; }

        for dj in 1..lgth2 + 1 {
            let gf2_j = gf2[dj];
            let gf2_jm1 = gf2[dj - 1];
            let diag = prev[dj - 1];

            let mut wm = diag;

            // mi (row-running gap-skip tracker)
            // C line 1328: g = mi + fgcp2[col_jm1] * gf1_i
            let g = fmadd(fgcp2[dj - 1], gf1_i, mi);
            if g > wm { wm = g; }
            // C line 1338: g = previousw[dj-1] + ogcp2[col_j] * gf1_im1
            let g = fmadd(ogcp2[dj], gf1_im1, diag);
            if g >= mi {
                mi = g;
                mpi = (dj - 1) as i64;
            }
            // C MSalignmm has `USE_PENALTY_EX = 0` (MSalignmm.c:7), so
            // the `mi += fpenalty_ex;` increment that `A__align`
            // applies is OMITTED here.

            // m[dj] (column-running gap-skip tracker)
            // C line 1349: g = m[dj] + fgcp1[row_im1] * gf2_j
            let g = fmadd(fgcp1_im1, gf2_j, m[dj]);
            if g > wm { wm = g; }
            // C line 1361: g = previousw[dj-1] + ogcp1[row_i] * gf2_jm1
            let g = fmadd(ogcp1_i, gf2_jm1, diag);
            if g >= m[dj] {
                m[dj] = g;
                mp[dj] = (di - 1) as i64;
            }
            // `USE_PENALTY_EX = 0`: skip `m[dj] += fpenalty_ex;` too.

            cur[dj] += wm;

            // At the midpoint, capture state for the join with backward.
            if capture {
                jumpbackj[dj] = mp[dj];
                jumpbacki[dj] = mpi;
                midw[dj] = cur[dj];
                midm[dj] = m[dj];
                midn[dj] = mi;
            }
        }
    }

    ForwardState { midw, midm, midn, jumpbackj, jumpbacki, currentw, previousw, m, mp }
}

/// Backward DP from row `ien` (= `ist + lgth1 - 1`) down to `imid - 1`.
/// Adds the partial scores to `fwd.midw / midm / midn` at the
/// midpoint, then finds the optimal split column. Mirrors C
/// `MSalignmm_rec:1442-1785`.
///
/// `tailgp` has no direct role here: the tail boundary enters through
/// `gf1[ien + 1]` / `gf2[jen + 1]`.
fn backward_dp(
    ctx: &mut MsCtx<'_>,
    ist: usize, ien: usize,
    jst: usize, jen: usize,
    imid: usize,
    fwd: ForwardState,
) -> BackwardState {
    let lgth1 = ien - ist + 1;
    let lgth2 = jen - jst + 1;
    debug_assert!(ien < ctx.prof1.length && jen < ctx.prof2.length);
    debug_assert!(lgth1 >= 2 && lgth2 >= 2);

    let tail1 = ctx.gf1[ien + 1];
    let tail2 = ctx.gf2[jen + 1];

    // initverticalw[di] = match_calc(prof2[jen], prof1[ist..=ien]) +
    //   fgcp1[ien] * gf2_lgth2 + ogcp1[ist+di+1] * gf2_lgth2m1
    let mut initverticalw = vec![0.0f64; lgth1 + 1];
    ctx.match_col(jen, ist, lgth1, &mut initverticalw[..lgth1]);
    {
        let gf2_tail = tail2;
        let gf2_lgth2m1 = ctx.gf2[jen];
        let fgcp1_end = ctx.fgcp1[ien];
        let ogcp1 = &ctx.ogcp1[ist + 1..=ien];
        for (di, iv) in initverticalw[..lgth1 - 1].iter_mut().enumerate() {
            *iv = fmadd(fgcp1_end, gf2_tail, fmadd(ogcp1[di], gf2_lgth2m1, *iv));
        }
    }

    // currentw[dj] = match_calc(prof1[ien], prof2[jst..=jen]) +
    //   fgcp2[jen] * gf1_lgth1 + ogcp2[jst+dj+1] * gf1_lgth1m1
    let ForwardState {
        mut midw, mut midm, mut midn, jumpbackj, jumpbacki,
        mut currentw, mut previousw, mut m, mut mp,
    } = fwd;
    currentw.fill(0.0);
    previousw.fill(0.0);
    ctx.match_row(ien, jst, lgth2, &mut currentw[..lgth2]);
    {
        let gf1_tail = tail1;
        let gf1_lgth1m1 = ctx.gf1[ien];
        let fgcp2_end = ctx.fgcp2[jen];
        let ogcp2 = &ctx.ogcp2[jst + 1..=jen];
        for (dj, cw) in currentw[..lgth2 - 1].iter_mut().enumerate() {
            *cw = fmadd(fgcp2_end, gf1_tail, fmadd(ogcp2[dj], gf1_lgth1m1, *cw));
        }
    }

    // m[dj] = currentw[dj+1] + fgcp1[ien-1] * gf2_{jen-dj+1}
    {
        let fgcp1_prev = ctx.fgcp1[ien - 1];
        let gf2 = &ctx.gf2[jst + 1..=jen + 1];
        for dj in (0..lgth2).rev() {
            m[dj] = fmadd(fgcp1_prev, gf2[dj], currentw[dj + 1]);
            mp[dj] = (lgth1 as i64) - 1;
        }
        m[lgth2] = 0.0;
        mp[lgth2] = 0;
    }

    // Storage for split-point detection.
    // C initialises `jumpi = 0` (`MSalignmm.c:1263`) and only overwrites
    // it at `i == imid - 1`. The backward loop then keeps going until
    // `i == jumpi` (`:1714`), so an m-split (`jumpi = jumpbackj[jmid]`,
    // which can sit well above the midpoint) records `jumpforwi` /
    // `jumpforwj` on that later row (`:1597`) and folds those rows into
    // `firstmp`. Stopping at the midpoint uses the wrong row's pointers
    // and can place the bottom half past one profile, dropping the
    // other profile's remaining columns.
    let mut jmid: usize = 0;
    let mut jumpi: usize = 0;
    let mut jumpj: usize = 0;
    let mut stop_i: i64 = 0;
    let mut split_ready = false;
    let mut jumpforwi = vec![0i64; lgth2 + 1];
    let mut jumpforwj = vec![0i64; lgth2 + 1];
    let mut firstm: f64 = f64::NEG_INFINITY;
    let mut firstmp: i64 = lgth1 as i64;

    let mut di_signed: i64 = (lgth1 as i64) - 2;
    while di_signed >= 0 {
        let di = di_signed as usize;
        std::mem::swap(&mut previousw, &mut currentw);
        previousw[lgth2 - 1] = initverticalw[di + 1];

        let row_i = ist + di;
        let row_ip1 = ist + di + 1;
        ctx.match_row(row_i, jst, lgth2, &mut currentw[..lgth2]);
        currentw[lgth2 - 1] = initverticalw[di];

        // Row-invariant gap terms.
        let gf1_i = ctx.gf1[row_i];
        let gf1_ip1 = ctx.gf1[row_ip1];
        let ogcp1_ip1 = ctx.ogcp1[row_ip1];
        let fgcp1_i = ctx.fgcp1[row_i];
        // Column windows, indexed by the local `dj`.
        let ogcp2 = &ctx.ogcp2[jst..=jen + 1];
        let fgcp2 = &ctx.fgcp2[jst..=jen];
        let gf2 = &ctx.gf2[jst..=jen + 1];
        let prev = &previousw[..=lgth2];
        let cur = &mut currentw[..=lgth2];
        let m = &mut m[..=lgth2];
        let mp = &mut mp[..=lgth2];

        let mut mi = fmadd(fgcp2[lgth2 - 2], gf1_ip1, prev[lgth2 - 1]);
        let mut mpi: i64 = (lgth2 - 1) as i64;

        let record_forw = di as i64 == stop_i || di + 1 == imid;
        let at_mid = di == imid;
        let at_mid_m1 = di + 1 == imid;

        let mut dj_signed: i64 = (lgth2 as i64) - 2;
        while dj_signed >= 0 {
            let dj = dj_signed as usize;
            let gf2_j = gf2[dj];
            let gf2_jp1 = gf2[dj + 1];

            // Diagonal default (C lines 1548-1550). In the C
            // MSalignmm backward inner loop, `*prept = previousw[j+1]`
            // because `prept` starts at `previousw + lgth2 - 1` and
            // decrements (so at iter j=lgth2-2, *prept = previousw[lgth2-1]
            // = previousw[(lgth2-2) + 1]). All four `*prept` references
            // below therefore use `previousw[dj + 1]`, NOT `previousw[dj]`.
            let pprev = prev[dj + 1];
            let mut wm = pprev;
            let mut ijpi: i64 = (di + 1) as i64;
            let mut ijpj: i64 = (dj + 1) as i64;

            // C line 1552 (mi candidate):
            //   g = mi + ogcp2[col_jp1] * gf1_i
            //   if g > wm: wm = g, ijpj = mpi, ijpi = i+1
            let g = fmadd(ogcp2[dj + 1], gf1_i, mi);
            if g > wm {
                wm = g;
                ijpj = mpi;
                ijpi = (di + 1) as i64;
            }

            // mi update (C line 1561). g = *prept + fgcp2[j] * gf1[i+1]
            //   = previousw[dj+1] + fgcp2[col_j] * gf1_ip1.
            let g = fmadd(fgcp2[dj], gf1_ip1, pprev);
            if g >= mi {
                mi = g;
                mpi = (dj + 1) as i64;
            }
            // C MSalignmm has `USE_PENALTY_EX = 0`: skip `mi += fpenalty_ex`.

            // C line 1575 (mj candidate):
            //   g = m[dj] + ogcp1[row_ip1] * gf2_j
            //   if g > wm: wm = g, ijpi = mp[dj], ijpj = j+1
            let g = fmadd(ogcp1_ip1, gf2_j, m[dj]);
            if g > wm {
                wm = g;
                ijpi = mp[dj];
                ijpj = (dj + 1) as i64;
            }

            // mj update (C line 1585). g = *prept + fgcp1[i] * gf2_{j+1}
            //   = previousw[dj+1] + fgcp1[row_i] * gf2_jp1.
            let g = fmadd(fgcp1_i, gf2_jp1, pprev);
            if g >= m[dj] {
                m[dj] = g;
                mp[dj] = (di + 1) as i64;
            }
            // `USE_PENALTY_EX = 0`: skip `m[dj] += fpenalty_ex`.

            // C `MSalignmm.c:1597`: record when `i == jumpi || i == imid-1`.
            // `stop_i` is still the value from before this row (0 until
            // the midpoint decision, then the m-split row), matching C
            // where `jumpi` is assigned after the inner j loop.
            if record_forw {
                jumpforwi[dj] = ijpi;
                jumpforwj[dj] = ijpj;
            }
            // Accumulate midw/midm at row imid; midn at row imid-1.
            // C `MSalignmm.c:1610-1612`:
            //   midw[j] += wm;        // NOTE: j, not j+1
            //   midm[j+1] += *mjpt;
            // The midw index off-by-one (j vs j+1) was a port bug
            // that surfaced as a 1-column drift on asymmetric inputs.
            if at_mid {
                midw[dj] += wm;
                midm[dj + 1] += m[dj];
            }
            if at_mid_m1 {
                midn[dj] += mi;
            }

            cur[dj] += wm;

            dj_signed -= 1;
        }
        // C line 1628: track firstm = max over rows of previousw[0] + fgcp1[i].
        let g_first = fgcp1_i + previousw[0];
        if firstm < g_first {
            firstm = g_first;
            firstmp = (di + 1) as i64;
        }
        // C line 1637: `if( i == imid ) midm[j+1] += firstm;` — at
        // this point the inner-loop j has decremented past 0 to -1,
        // so `j + 1` is 0. C touches `midm[0]`.
        if at_mid {
            midm[0] += firstm;
        }

        // At i == imid - 1, decide jmid + (jumpi, jumpj). Do not break:
        // C breaks only when `i == jumpi` (`MSalignmm.c:1714`), which for
        // an m-split is a later row.
        if !split_ready && di + 1 == imid {
            let mut maxwm = midw.get(1).copied().unwrap_or(f64::NEG_INFINITY);
            jmid = 0;
            for j in 2..lgth2.saturating_sub(1) {
                if let Some(&w) = midw.get(j) {
                    if w > maxwm { maxwm = w; jmid = j; }
                }
            }
            for j in 0..=lgth2 {
                if let Some(&w) = midm.get(j) {
                    if w > maxwm { maxwm = w; jmid = j; }
                }
            }

            // Reconcile: which of {midw, midm, midn} produced max?
            let wmw = midw.get(jmid).copied().unwrap_or(f64::NEG_INFINITY);
            jumpi = imid.saturating_sub(1);
            jumpj = jmid.saturating_sub(1);
            let mut wmsel = wmw;
            if jmid > 0 {
                let nval = midn.get(jmid.saturating_sub(1)).copied().unwrap_or(f64::NEG_INFINITY);
                if nval > wmsel {
                    jumpi = imid.saturating_sub(1);
                    let jb = jumpbacki.get(jmid).copied().unwrap_or(0);
                    jumpj = if jb < 0 { 0 } else { jb as usize };
                    wmsel = nval;
                }
            }
            let mval = midm.get(jmid).copied().unwrap_or(f64::NEG_INFINITY);
            if mval > wmsel {
                let jb = jumpbackj.get(jmid).copied().unwrap_or(0);
                jumpi = if jb < 0 { 0 } else { jb as usize };
                jumpj = jmid.saturating_sub(1);
            }
            stop_i = jumpi as i64;
            split_ready = true;
        }

        // C lines 1714-1770, on the row `i == jumpi`. `firstmp` above
        // already includes this row.
        if split_ready && di as i64 == stop_i {
            let mut imid_eff = imid;
            if jmid == 0 {
                // C: `if( imid < firstmp-1 ) { jumpi = firstmp; imid = firstmp+1; }`
                if (imid as i64) < firstmp - 1 && firstmp >= 0 {
                    jumpi = firstmp as usize;
                    imid_eff = firstmp as usize + 1;
                }
                jmid = 1;
                jumpj = 0;
            } else if jmid >= lgth2 {
                jumpi = imid.saturating_sub(1);
                jmid = lgth2;
                jumpj = lgth2.saturating_sub(1);
            } else {
                // `imid = jumpforwi[jumpj]; jmid = jumpforwj[jumpj]`
                // from the row we just finished, which is `jumpi` and
                // not necessarily `imid - 1`.
                let new_imid = jumpforwi.get(jumpj).copied().unwrap_or(imid as i64);
                let new_jmid = jumpforwj.get(jumpj).copied().unwrap_or(jmid as i64);
                if new_imid >= 0 && new_jmid >= 0 {
                    imid_eff = new_imid as usize;
                    jmid = new_jmid as usize;
                }
                if imid_eff == jumpi {
                    jumpi = imid_eff.saturating_sub(1);
                }
            }
            return BackwardState {
                split_point: (jumpi, jumpj, jmid),
                imid_override: Some(imid_eff),
            };
        }

        di_signed -= 1;
    }

    BackwardState {
        split_point: (jumpi, jumpj, jmid),
        imid_override: Some(imid),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::profile::Profile;
    use crate::dp::GapModel;

    fn simple_setup() -> (Vec<Vec<f64>>, [u8; 256]) {
        let mut mtx = vec![vec![-100.0f64; 5]; 5];
        for i in 0..4 { mtx[i][i] = 100.0; }
        let mut map = [4u8; 256];
        map[b'A' as usize] = 0;
        map[b'C' as usize] = 1;
        map[b'G' as usize] = 2;
        map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;
        (mtx, map)
    }

    /// On a SMALL input (< DPTANNI), msalignmm must produce output
    /// byte-identical to `profile_align_imp_with_boundary` — because
    /// the base case just delegates to it.
    #[test]
    fn base_case_matches_full_dp() {
        let (mtx, map) = simple_setup();
        let s1: &[u8] = b"ACGTACGTAC";
        let s2: &[u8] = b"ACGTACGTAC";
        let p1 = Profile::from_aligned(&[s1], &[1.0], &map, 5);
        let p2 = Profile::from_aligned(&[s2], &[1.0], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);

        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, true, true, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, true, true);

        assert!((aln_full.score - aln_ms.score).abs() < 1e-9,
            "scores differ: full={} ms={}", aln_full.score, aln_ms.score);
        assert_eq!(aln_full.operations, aln_ms.operations,
            "trace differs for short input");
    }

    /// Two 360-char identical-prefix sequences (length similar to
    /// the rhodopsin sample). Should align diagonally end-to-end
    /// with `head_gap=false, tail_gap=false`. If msalignmm gives a
    /// shorter trace, the recursion is broken at this length.
    #[test]
    fn freegap_360char_diagonal() {
        let (mtx, map) = simple_setup();
        let s: Vec<u8> = (0..360).map(|i| b"ACGT"[(i * 3 + 1) % 4]).collect();
        let p1 = Profile::from_aligned(&[s.as_slice()], &[1.0], &map, 5);
        let p2 = Profile::from_aligned(&[s.as_slice()], &[1.0], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);
        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, false, false, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, false, false);
        assert_eq!(aln_full.operations.len(), aln_ms.operations.len(),
            "360-char free-gap: trace length differs: full={} ms={}",
            aln_full.operations.len(), aln_ms.operations.len());
        assert!((aln_full.score - aln_ms.score).abs() < 1e-6);
    }

    /// Asymmetric multi-seq Hirschberg (100 chars × 105-char 2-seq
    /// profile). Mirrors a typical progressive merge: a single new
    /// sequence joined to a small existing group.
    #[test]
    fn freegap_100_vs_105_multiseq() {
        let (mtx, map) = simple_setup();
        let s1: Vec<u8> = (0..100).map(|i| b"ACGT"[(i * 3 + 1) % 4]).collect();
        let s2_g1: Vec<u8> = (0..105).map(|i| b"ACGT"[(i * 7 + 5) % 4]).collect();
        let s2_g2: Vec<u8> = (0..105).map(|i| b"ACGT"[(i * 11 + 3) % 4]).collect();
        let p1 = Profile::from_aligned(&[s1.as_slice()], &[1.0], &map, 5);
        let p2 = Profile::from_aligned(&[s2_g1.as_slice(), s2_g2.as_slice()], &[0.5, 0.5], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);
        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, true, true, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, true, true);
        assert_eq!(aln_full.operations.len(), aln_ms.operations.len(),
            "multi-seq asymmetric: trace length differs");
        assert!((aln_full.score - aln_ms.score).abs() < 1e-6);
    }

    /// MULTI-sequence profile case (what the progressive merge
    /// actually feeds msalignmm). Verifies msalignmm handles
    /// profiles built from 2+ aligned sequences (with internal
    /// gap columns and non-trivial nongap_freq distributions).
    #[test]
    fn multiseq_profile_matches_full_dp() {
        let (mtx, map) = simple_setup();
        // Pre-aligned group of 3 sequences with internal gaps.
        let group1: Vec<&[u8]> = vec![
            b"AC-GTACGTAC-GTAC",
            b"ACGGTAC-TACGGT-C",
            b"AC-GTACGTACGGTAC",
        ];
        let group2: Vec<&[u8]> = vec![
            b"ACGTACGTAC-GTAC",
            b"ACGGTACGT-CGTAC",
        ];
        let p1 = Profile::from_aligned(&group1, &[1.0/3.0; 3], &map, 5);
        let p2 = Profile::from_aligned(&group2, &[0.5, 0.5], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);
        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, false, false, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, false, false);
        assert!((aln_full.score - aln_ms.score).abs() < 1e-6,
            "multi-seq profile: scores differ: full={} ms={}",
            aln_full.score, aln_ms.score);
        assert_eq!(aln_full.operations, aln_ms.operations,
            "multi-seq profile: trace differs");
    }

    /// `head_gap=false, tail_gap=false` case (FFT-NS-2 / NW-NS-2 path
    /// in the engine). Verifies msalignmm matches the full DP for
    /// the free-terminal-gap settings the progressive merge uses.
    #[test]
    fn freegap_matches_full_dp() {
        let (mtx, map) = simple_setup();
        let s: Vec<u8> = (0..200).map(|i| b"ACGT"[(i * 3 + 1) % 4]).collect();
        let s2: Vec<u8> = (0..200).map(|i| b"ACGT"[(i * 7 + 2) % 4]).collect();
        let p1 = Profile::from_aligned(&[s.as_slice()], &[1.0], &map, 5);
        let p2 = Profile::from_aligned(&[s2.as_slice()], &[1.0], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);
        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, false, false, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, false, false);
        assert!((aln_full.score - aln_ms.score).abs() < 1e-6,
            "head_gap=false: scores differ: full={} ms={}", aln_full.score, aln_ms.score);
        assert_eq!(aln_full.operations.len(), aln_ms.operations.len(),
            "head_gap=false: trace length differs: full={} ms={}",
            aln_full.operations.len(), aln_ms.operations.len());
    }

    /// Recursive case with a non-trivial gap-required alignment.
    /// We embed a UNIQUE motif on either side of the gap region so
    /// the optimal alignment is uniquely determined (no tie-breaks).
    /// If the Hirschberg split logic is wrong, the score will drop.
    #[test]
    fn recursive_case_with_gap_matches_full_dp() {
        // Use the full BLOSUM-like 5-letter alphabet, with distinct
        // motifs anchoring the gap region.
        let mut mtx = vec![vec![-100.0f64; 5]; 5];
        for i in 0..4 { mtx[i][i] = 100.0; }
        let mut map = [4u8; 256];
        map[b'A' as usize] = 0;
        map[b'C' as usize] = 1;
        map[b'G' as usize] = 2;
        map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;

        // Left half: random-ish, right half: distinct motif — gives a
        // unique optimal alignment when 20-residue insert is added.
        let left: Vec<u8> = (0..60).map(|i| b"ACGT"[(i * 3 + 1) % 4]).collect();
        let right: Vec<u8> = (0..60).map(|i| b"TGCA"[(i * 7 + 5) % 4]).collect();
        let insert: Vec<u8> = (0..20).map(|i| b"AAGG"[(i + 2) % 4]).collect();

        let s_short: Vec<u8> = left.iter().chain(right.iter()).copied().collect();
        let s_long: Vec<u8> = left.iter()
            .chain(insert.iter())
            .chain(right.iter())
            .copied().collect();
        let p1 = Profile::from_aligned(&[s_short.as_slice()], &[1.0], &map, 5);
        let p2 = Profile::from_aligned(&[s_long.as_slice()], &[1.0], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);

        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, true, true, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, true, true);

        assert!((aln_full.score - aln_ms.score).abs() < 1e-6,
            "scores differ: full={} ms={}",
            aln_full.score, aln_ms.score);
        // Trace must match too — input has a unique optimal alignment
        // (distinct anchor motifs on either side of the gap).
        assert_eq!(aln_full.operations, aln_ms.operations,
            "trace differs for gap-required input with unique optimum");
    }

    /// Recursive case: an input long enough to trigger Hirschberg
    /// (lgth1 >= DPTANNI=100 AND lgth2 >= DPTANNI). The trace must
    /// still match the full DP, because Hirschberg gives the same
    /// optimal alignment.
    #[test]
    fn recursive_case_matches_full_dp() {
        let (mtx, map) = simple_setup();
        // Identical 200-nt sequences — the optimal alignment is
        // diagonal everywhere; Hirschberg must reproduce it exactly.
        let s: Vec<u8> = (0..200).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let p1 = Profile::from_aligned(&[s.as_slice()], &[1.0], &map, 5);
        let p2 = Profile::from_aligned(&[s.as_slice()], &[1.0], &map, 5);
        let gap = GapModel::new(-200.0, -10.0);

        let aln_full = profile_align_imp_with_boundary(
            &p1, &p2, &mtx, &gap, true, true, None, false,
            BoundaryFreqs::default(),
        );
        let aln_ms = msalignmm(&p1, &p2, &mtx, &gap, true, true);

        assert!((aln_full.score - aln_ms.score).abs() < 1.0,
            "scores differ: full={} ms={}", aln_full.score, aln_ms.score);
        // The trace MUST be byte-identical for identical inputs.
        assert_eq!(aln_full.operations.len(), aln_ms.operations.len(),
            "trace length differs: full={} ms={}",
            aln_full.operations.len(), aln_ms.operations.len());
        assert_eq!(aln_full.operations, aln_ms.operations,
            "trace differs for identical 200-nt input");
    }

    /// Gappy multi-sequence profiles long enough to recurse. The
    /// Hirschberg split must keep every input column: a `jmid == lgth2`
    /// or m-split used to return a trace narrower than the longer
    /// profile.
    #[test]
    fn gappy_recursive_profiles_keep_every_column() {
        let (mtx, map) = simple_setup();
        let gap = GapModel::new(-1200.0, -60.0);
        let pattern = |n: usize, salt: usize, gap_every: usize| -> Vec<u8> {
            (0..n).map(|i| {
                if gap_every > 0 && (i + salt) % gap_every == 0 { b'-' } else { b"ACGT"[(i * salt + 3) % 4] }
            }).collect()
        };
        for (n1, n2, g1, g2) in [(180, 220, 7, 11), (250, 140, 5, 9), (300, 300, 4, 13)] {
            let a = pattern(n1, 1, g1);
            let b = pattern(n1, 2, g1);
            let c = pattern(n2, 3, g2);
            let d = pattern(n2, 5, g2);
            let p1 = Profile::from_aligned(&[a.as_slice(), b.as_slice()], &[0.5, 0.5], &map, 5);
            let p2 = Profile::from_aligned(&[c.as_slice(), d.as_slice()], &[0.5, 0.5], &map, 5);
            let aln = msalignmm(&p1, &p2, &mtx, &gap, true, true);
            let kept1 = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match | AlignOp::Delete)).count();
            let kept2 = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match | AlignOp::Insert)).count();
            assert_eq!(kept1, p1.length, "dropped profile-1 columns at {n1}x{n2}");
            assert_eq!(kept2, p2.length, "dropped profile-2 columns at {n1}x{n2}");
            assert!(aln.operations.len() >= p1.length.max(p2.length));
        }
    }
}
