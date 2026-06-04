/*
 * Instrumented copy of MSalignmm_rec's first-level forward + backward
 * DP, exposing midw[], midm[], midn[], jumpbacki/j[], jumpforwi/j[],
 * and the chosen (jmid, jumpi, jumpj) at the top-level split decision.
 *
 * Used by `crates/mafft-core/tests/cross_validate_msalign.rs` to
 * pinpoint where our Rust `msalignmm` Hirschberg DP diverges from C
 * `MSalignmm_rec`. Re-implements just the top-level recursion (no
 * tail recursion) — every cell-update mirrors C MSalignmm.c
 * exactly.
 *
 * Symbols are prefixed `rs_msi_` so they don't clash with the
 * upstream library.
 */

#include "mltaln.h"
#include "dp.h"

#define MEMSAVE 1
#define USE_PENALTY_EX 0
#define FASTMATCHCALC 1

/* Match-score helper: same as MSalignmm.c::match_calc_alphabet_seq with FASTMATCHCALC=1. */
static void rs_msi_match_calc(
    double **n_dynamicmtx, double *match,
    double **cpmx1, double **cpmx2,
    int i1, int start2, int lgth2,
    double **doublework, int **intwork, int initialize)
{
    int j, l, p;
    double **cpmxpd = doublework;
    int **cpmxpdn = intwork;
    double *matchpt, *cpmxpdpt, **cpmxpdptpt;
    int *cpmxpdnpt, **cpmxpdnptpt;
    double *scarr;
    scarr = calloc(nalphabets, sizeof(double));

    if (initialize) {
        int count = 0;
        for (j = 0, p = start2; j < lgth2; j++, p++) {
            count = 0;
            for (l = 0; l < nalphabets; l++) {
                if (cpmx2[l][p]) {
                    cpmxpd[j][count] = cpmx2[l][p];
                    cpmxpdn[j][count] = l;
                    count++;
                }
            }
            cpmxpdn[j][count] = -1;
        }
    }

    for (l = 0; l < nalphabets; l++) {
        scarr[l] = 0.0;
        for (j = 0; j < nalphabets; j++)
            scarr[l] += n_dynamicmtx[j][l] * cpmx1[j][i1];
    }
    matchpt = match;
    cpmxpdnptpt = cpmxpdn;
    cpmxpdptpt = cpmxpd;
    while (lgth2--) {
        *matchpt = 0.0;
        cpmxpdnpt = *cpmxpdnptpt++;
        cpmxpdpt = *cpmxpdptpt++;
        while (*cpmxpdnpt > -1)
            *matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
        matchpt++;
    }
    free(scarr);
}

/*
 * Run the first-level forward + backward DP at the top of
 * MSalignmm_rec (no recursion), capture midw/midm/midn,
 * jumpback*, jumpforw*, and the split decision.
 *
 * Sequences: single-sequence groups (icyc = jcyc = 1).
 * Caller allocates all output arrays to size lgth2 + 2.
 */
void rs_msalignmm_capture_top(
    double **n_dynamicmtx,
    char *seq1, char *seq2,
    int lgth1, int lgth2,
    int headgp, int tailgp,
    /* outputs */
    int *out_imid, int *out_jmid, int *out_jumpi, int *out_jumpj,
    double *out_midw,
    double *out_midm,
    double *out_midn,
    int *out_jumpbacki,
    int *out_jumpbackj,
    int *out_jumpforwi,
    int *out_jumpforwj)
{
    int i, j;
    int ll1 = lgth1 + 100, ll2 = lgth2 + 100;
    int imid = lgth1 / 2;
    double wm;

    /* Build single-sequence inputs. */
    char *seqs1[1] = {seq1};
    char *seqs2[1] = {seq2};
    double eff1[1] = {1.0};
    double eff2[1] = {1.0};

    /* cpmx matrices: [alphabet][position]. */
    double **cpmx1 = AllocateFloatMtx(nalphabets, ll1 + 2);
    double **cpmx2 = AllocateFloatMtx(nalphabets, ll2 + 2);
    cpmx_calc_new(seqs1, cpmx1, eff1, lgth1, 1);
    cpmx_calc_new(seqs2, cpmx2, eff2, lgth2, 1);

    /* Gap counts. */
    double *ogcp1opt = AllocateFloatVec(ll1 + 2);
    double *ogcp2opt = AllocateFloatVec(ll2 + 2);
    double *fgcp1opt = AllocateFloatVec(ll1 + 2);
    double *fgcp2opt = AllocateFloatVec(ll2 + 2);
    st_OpeningGapCount(ogcp1opt, 1, seqs1, eff1, lgth1);
    st_FinalGapCount(fgcp1opt, 1, seqs1, eff1, lgth1);
    st_OpeningGapCount(ogcp2opt, 1, seqs2, eff2, lgth2);
    st_FinalGapCount(fgcp2opt, 1, seqs2, eff2, lgth2);

    double *gapfreq1f = AllocateFloatVec(ll1 + 2);
    double *gapfreq2f = AllocateFloatVec(ll2 + 2);
    gapcountf(gapfreq1f, seqs1, 1, eff1, lgth1);
    gapcountf(gapfreq2f, seqs2, 1, eff2, lgth2);
    for (i = 0; i < lgth1 + 1; i++) gapfreq1f[i] = 1.0 - gapfreq1f[i];
    for (i = 0; i < lgth2 + 1; i++) gapfreq2f[i] = 1.0 - gapfreq2f[i];
    double headgapfreq1 = 1.0;
    double headgapfreq2 = 1.0;

    /* Final ogcp/fgcp arrays (scaled). */
    double fpenalty = (double)penalty;
    double *ogcp1 = AllocateFloatVec(ll1 + 2);
    double *ogcp2 = AllocateFloatVec(ll2 + 2);
    double *fgcp1 = AllocateFloatVec(ll1 + 2);
    double *fgcp2 = AllocateFloatVec(ll2 + 2);
    for (i = 0; i < lgth1; i++) {
        ogcp1[i] = 0.5 * (1.0 - ogcp1opt[i]) * fpenalty * gapfreq1f[i];
        fgcp1[i] = 0.5 * (1.0 - fgcp1opt[i]) * fpenalty * gapfreq1f[i];
    }
    for (i = 0; i < lgth2; i++) {
        ogcp2[i] = 0.5 * (1.0 - ogcp2opt[i]) * fpenalty * gapfreq2f[i];
        fgcp2[i] = 0.5 * (1.0 - fgcp2opt[i]) * fpenalty * gapfreq2f[i];
    }

    /* DP scratch arrays (sizes match C MSalignmm_rec). */
    double *w1 = AllocateFloatVec(ll2 + 2);
    double *w2 = AllocateFloatVec(ll2 + 2);
    double *midw = AllocateFloatVec(ll2 + 2);
    double *midn = AllocateFloatVec(ll2 + 2);
    double *midm = AllocateFloatVec(ll2 + 2);
    int *jumpbacki = AllocateIntVec(ll2 + 2);
    int *jumpbackj = AllocateIntVec(ll2 + 2);
    int *jumpforwi = AllocateIntVec(ll2 + 2);
    int *jumpforwj = AllocateIntVec(ll2 + 2);

    double *initverticalw = AllocateFloatVec(ll1 + 2);
    double *lastverticalw = AllocateFloatVec(ll1 + 2);

    double *m = AllocateFloatVec(ll2 + 2);
    int *mp = AllocateIntVec(ll2 + 2);

    double **doublework = AllocateFloatMtx(ll1 > ll2 ? ll1 + 2 : ll2 + 2, nalphabets);
    int **intwork = AllocateIntMtx(ll1 > ll2 ? ll1 + 2 : ll2 + 2, nalphabets);

    double *currentw = w1;
    double *previousw = w2;

    /* --- Forward DP --- */
    rs_msi_match_calc(n_dynamicmtx, initverticalw, cpmx2, cpmx1, 0, 0, lgth1, doublework, intwork, 1);
    rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, 0, 0, lgth2, doublework, intwork, 1);

    for (i = 1; i < lgth1 + 1; i++)
        initverticalw[i] += (ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2f[0]);
    for (j = 1; j < lgth2 + 1; j++)
        currentw[j] += (ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1f[0]);
    (void)headgp; /* always true for our test */
    (void)tailgp;

    for (j = 1; j < lgth2 + 1; ++j) {
        m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2f[j - 1];
        mp[j] = 0;
    }
    lastverticalw[0] = currentw[lgth2 - 1];

    for (i = 1; i <= imid; i++) {
        double *wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;
        previousw[0] = initverticalw[i - 1];

        rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, i, 0, lgth2, doublework, intwork, 0);
        currentw[0] = initverticalw[i];
        m[0] = ogcp1[i];
        if (i == imid) midm[0] = m[0];

        double mi = previousw[0] + ogcp2[1] * gapfreq1f[i - 1];
        int mpi = 0;

        double *mjpt = m + 1;
        double *prept = previousw;
        double *curpt = currentw + 1;
        int *mpjpt = mp + 1;

        for (j = 1; j < lgth2 + 1; j++) {
            wm = *prept;
            double g = mi + fgcp2[j - 1] * gapfreq1f[i];
            if (g > wm) wm = g;
            g = *prept + ogcp2[j] * gapfreq1f[i - 1];
            if (g >= mi) { mi = g; mpi = j - 1; }
            g = *mjpt + fgcp1[i - 1] * gapfreq2f[j];
            if (g > wm) wm = g;
            g = *prept + ogcp1[i] * gapfreq2f[j - 1];
            if (g >= *mjpt) { *mjpt = g; *mpjpt = i - 1; }
            *curpt += wm;

            if (i == imid) {
                jumpbackj[j] = *mpjpt;
                jumpbacki[j] = mpi;
                midw[j] = *curpt;
                midm[j] = *mjpt;
                midn[j] = mi;
            }
            mjpt++; prept++; mpjpt++; curpt++;
        }
        lastverticalw[i] = currentw[lgth2 - 1];
    }

    /* --- Backward DP --- */
    rs_msi_match_calc(n_dynamicmtx, initverticalw, cpmx2, cpmx1, lgth2 - 1, 0, lgth1, doublework, intwork, 1);
    rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, lgth1 - 1, 0, lgth2, doublework, intwork, 1);

    for (i = 0; i < lgth1 - 1; i++)
        initverticalw[i] += (fgcp1[lgth1 - 1] * gapfreq2f[lgth2] + ogcp1[i + 1] * gapfreq2f[lgth2 - 1]);
    for (j = 0; j < lgth2 - 1; j++)
        currentw[j] += (fgcp2[lgth2 - 1] * gapfreq1f[lgth1] + ogcp2[j + 1] * gapfreq1f[lgth1 - 1]);

    for (j = lgth2 - 1; j > -1; --j) {
        m[j] = currentw[j + 1] + fgcp1[lgth1 - 2] * gapfreq2f[j + 1];
        mp[j] = lgth1 - 1;
    }

    double firstm = -9999999.9;
    int firstmp = lgth1;
    int jumpi = 0, jumpj = 0, jmid = 0;
    double maxwm = 0.0;

    for (i = lgth1 - 2; i > -1; i--) {
        double *wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;
        previousw[lgth2 - 1] = initverticalw[i + 1];

        rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, i, 0, lgth2, doublework, intwork, 0);
        currentw[lgth2 - 1] = initverticalw[i];

        double mi = previousw[lgth2 - 1] + fgcp2[lgth2 - 2] * gapfreq1f[i + 1];
        int mpi = lgth2 - 1;

        double *mjpt = m + lgth2 - 2;
        double *prept = previousw + lgth2 - 1;
        double *curpt = currentw + lgth2 - 2;
        int *mpjpt = mp + lgth2 - 2;

        for (j = lgth2 - 2; j > -1; j--) {
            wm = *prept;
            int ijpi = i + 1;
            int ijpj = j + 1;

            double g = mi + ogcp2[j + 1] * gapfreq1f[i];
            if (g > wm) { wm = g; ijpj = mpi; ijpi = i + 1; }

            g = *prept + fgcp2[j] * gapfreq1f[i + 1];
            if (g >= mi) { mi = g; mpi = j + 1; }

            g = *mjpt + ogcp1[i + 1] * gapfreq2f[j];
            if (g > wm) { wm = g; ijpi = *mpjpt; ijpj = j + 1; }

            g = *prept + fgcp1[i] * gapfreq2f[j + 1];
            if (g >= *mjpt) { *mjpt = g; *mpjpt = i + 1; }

            if (i == jumpi || i == imid - 1) {
                jumpforwi[j] = ijpi;
                jumpforwj[j] = ijpj;
            }
            if (i == imid) {
                midw[j + 1] += wm;
                midm[j + 1] += *mjpt;
            }
            if (i == imid - 1) midn[j] += mi;

            *curpt += wm;
            mjpt--; prept--; mpjpt--; curpt--;
        }
        double g = *prept + fgcp1[i];
        if (firstm < g) { firstm = g; firstmp = i + 1; }
        if (i == imid) midm[j + 1] += firstm;

        if (i == imid - 1) {
            maxwm = midw[1];
            jmid = 0;
            for (j = 2; j < lgth2 - 1; j++) {
                wm = midw[j];
                if (wm > maxwm) { jmid = j; maxwm = wm; }
            }
            for (j = 0; j < lgth2 + 1; j++) {
                wm = midm[j];
                if (wm > maxwm) { jmid = j; maxwm = wm; }
            }
            wm = midw[jmid];
            jumpi = imid - 1;
            jumpj = jmid - 1;
            if (jmid > 0 && midn[jmid - 1] > wm) {
                jumpi = imid - 1;
                jumpj = jumpbacki[jmid];
                wm = midn[jmid - 1];
            }
            if (midm[jmid] > wm) {
                jumpi = jumpbackj[jmid];
                jumpj = jmid - 1;
                wm = midm[jmid];
            }
            break;
        }
    }

    /* Edge cases (C lines 1721-1770). */
    if (jmid == 0) {
        if (imid < firstmp - 1) {
            jumpi = firstmp;
            imid = firstmp + 1;
        }
        jumpj = 0;
        jmid = 1;
    } else if (jmid >= lgth2) {
        jumpi = imid - 1;
        jmid = lgth2;
        jumpj = lgth2 - 1;
    } else {
        imid = jumpforwi[jumpj];
        jmid = jumpforwj[jumpj];
        if (imid == jumpi) jumpi = imid - 1;
    }

    *out_imid = imid;
    *out_jmid = jmid;
    *out_jumpi = jumpi;
    *out_jumpj = jumpj;
    for (j = 0; j < lgth2 + 1; j++) {
        out_midw[j] = midw[j];
        out_midm[j] = midm[j];
        out_midn[j] = midn[j];
        out_jumpbacki[j] = jumpbacki[j];
        out_jumpbackj[j] = jumpbackj[j];
        out_jumpforwi[j] = jumpforwi[j];
        out_jumpforwj[j] = jumpforwj[j];
    }

    /* Free everything. */
    FreeFloatMtx(cpmx1);
    FreeFloatMtx(cpmx2);
    FreeFloatVec(ogcp1opt);
    FreeFloatVec(ogcp2opt);
    FreeFloatVec(fgcp1opt);
    FreeFloatVec(fgcp2opt);
    FreeFloatVec(gapfreq1f);
    FreeFloatVec(gapfreq2f);
    FreeFloatVec(ogcp1);
    FreeFloatVec(ogcp2);
    FreeFloatVec(fgcp1);
    FreeFloatVec(fgcp2);
    FreeFloatVec(w1);
    FreeFloatVec(w2);
    FreeFloatVec(midw);
    FreeFloatVec(midn);
    FreeFloatVec(midm);
    FreeIntVec(jumpbacki);
    FreeIntVec(jumpbackj);
    FreeIntVec(jumpforwi);
    FreeIntVec(jumpforwj);
    FreeFloatVec(initverticalw);
    FreeFloatVec(lastverticalw);
    FreeFloatVec(m);
    FreeIntVec(mp);
    FreeFloatMtx(doublework);
    FreeIntMtx(intwork);
}

/*
 * Run C MSalignmm_tanni for a SUB-REGION of the parent profile.
 * Sets up cpmx and gapinfo as the parent MSalignmm would for the
 * FULL input, then calls our re-implementation of MSalignmm_tanni
 * with the sub-region indices (ist, ien, jst, jen). Returns the
 * aligned strings via the out_seq1/out_seq2 buffers (caller must
 * allocate to at least `ien - ist + jen - jst + 100` bytes).
 *
 * For SUB-REGION testing: caller passes the FULL parent seq1/seq2,
 * lgth1/lgth2 (= parent lengths), and ist/ien/jst/jen identifying
 * the sub-region within those parent sequences.
 *
 * Use this to verify whether `profile_align_imp_with_boundary` on
 * the equivalent Rust sub-profile produces the same output as C
 * MSalignmm_tanni in-context (where gapinfo arrays use PARENT
 * positions, possibly differing from a standalone build).
 */
void rs_msalignmm_tanni_capture(
    double **n_dynamicmtx,
    char *seq1, char *seq2,
    int lgth1, int lgth2,
    int ist, int ien, int jst, int jen,
    int headgp, int tailgp,
    char *out_seq1, char *out_seq2, int *out_width)
{
    int i, j;
    int sub_lgth1 = ien - ist + 1;
    int sub_lgth2 = jen - jst + 1;
    int ll1 = lgth1 + 100, ll2 = lgth2 + 100;
    int sub_ll1 = sub_lgth1 + 100, sub_ll2 = sub_lgth2 + 100;

    char *seqs1[1] = {seq1};
    char *seqs2[1] = {seq2};
    double eff1[1] = {1.0};
    double eff2[1] = {1.0};

    /* Build PARENT cpmx / gap arrays. */
    double **cpmx1 = AllocateFloatMtx(nalphabets, ll1 + 2);
    double **cpmx2 = AllocateFloatMtx(nalphabets, ll2 + 2);
    cpmx_calc_new(seqs1, cpmx1, eff1, lgth1, 1);
    cpmx_calc_new(seqs2, cpmx2, eff2, lgth2, 1);

    double *ogcp1opt = AllocateFloatVec(ll1 + 2);
    double *ogcp2opt = AllocateFloatVec(ll2 + 2);
    double *fgcp1opt = AllocateFloatVec(ll1 + 2);
    double *fgcp2opt = AllocateFloatVec(ll2 + 2);
    st_OpeningGapCount(ogcp1opt, 1, seqs1, eff1, lgth1);
    st_FinalGapCount(fgcp1opt, 1, seqs1, eff1, lgth1);
    st_OpeningGapCount(ogcp2opt, 1, seqs2, eff2, lgth2);
    st_FinalGapCount(fgcp2opt, 1, seqs2, eff2, lgth2);

    double *gapfreq1f = AllocateFloatVec(ll1 + 2);
    double *gapfreq2f = AllocateFloatVec(ll2 + 2);
    gapcountf(gapfreq1f, seqs1, 1, eff1, lgth1);
    gapcountf(gapfreq2f, seqs2, 1, eff2, lgth2);
    for (i = 0; i < lgth1 + 1; i++) gapfreq1f[i] = 1.0 - gapfreq1f[i];
    for (i = 0; i < lgth2 + 1; i++) gapfreq2f[i] = 1.0 - gapfreq2f[i];
    double headgapfreq1 = 1.0;
    double headgapfreq2 = 1.0;

    double fpenalty = (double)penalty;
    double *ogcp1 = AllocateFloatVec(ll1 + 2);
    double *ogcp2 = AllocateFloatVec(ll2 + 2);
    double *fgcp1 = AllocateFloatVec(ll1 + 2);
    double *fgcp2 = AllocateFloatVec(ll2 + 2);
    for (i = 0; i < lgth1; i++) {
        ogcp1[i] = 0.5 * (1.0 - ogcp1opt[i]) * fpenalty * gapfreq1f[i];
        fgcp1[i] = 0.5 * (1.0 - fgcp1opt[i]) * fpenalty * gapfreq1f[i];
    }
    for (i = 0; i < lgth2; i++) {
        ogcp2[i] = 0.5 * (1.0 - ogcp2opt[i]) * fpenalty * gapfreq2f[i];
        fgcp2[i] = 0.5 * (1.0 - fgcp2opt[i]) * fpenalty * gapfreq2f[i];
    }

    /* Slice gap arrays per MSalignmm_tanni:709-712. */
    double *sub_ogcp1 = ogcp1 + ist;
    double *sub_fgcp1 = fgcp1 + ist;
    double *sub_ogcp2 = ogcp2 + jst;
    double *sub_fgcp2 = fgcp2 + jst;
    double *sub_gapfreq1f = gapfreq1f + ist;
    double *sub_gapfreq2f = gapfreq2f + jst;
    double sub_headgapfreq1 = (ist > 0) ? sub_gapfreq1f[-1] : headgapfreq1;
    double sub_headgapfreq2 = (jst > 0) ? sub_gapfreq2f[-1] : headgapfreq2;

    /* Sliced cpmx pointers: cpmx[a][parent_pos] is contiguous, so
     * sub access is just cpmx[a][ist + sub_pos]. We pass cpmx as
     * the FULL matrix and rely on match_calc starting at sub-indices. */

    /* Allocate DP scratch (sub-region sized). */
    double *w1 = AllocateFloatVec(sub_ll2 + 2);
    double *w2 = AllocateFloatVec(sub_ll2 + 2);
    double *initverticalw = AllocateFloatVec(sub_ll1 + 2);
    double *lastverticalw = AllocateFloatVec(sub_ll1 + 2);
    double *m = AllocateFloatVec(sub_ll2 + 2);
    int *mp = AllocateIntVec(sub_ll2 + 2);
    double **doublework = AllocateFloatMtx(sub_ll1 > sub_ll2 ? sub_ll1 + 2 : sub_ll2 + 2, nalphabets + 1);
    int **intwork = AllocateIntMtx(sub_ll1 > sub_ll2 ? sub_ll1 + 2 : sub_ll2 + 2, nalphabets + 1);
    int **intmtx = AllocateIntMtx(sub_ll1 + 1, sub_ll2 + 1);
    int **ijp = intmtx;

    double *currentw = w1;
    double *previousw = w2;

    /* DP — exact copy of MSalignmm_tanni:776-888 with sliced arrays. */
    rs_msi_match_calc(n_dynamicmtx, initverticalw, cpmx2, cpmx1, jst, ist, sub_lgth1, doublework, intwork, 1);
    rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist, jst, sub_lgth2, doublework, intwork, 1);

    if (headgp || ist != 0) {
        for (i = 1; i < sub_lgth1 + 1; i++)
            initverticalw[i] += (sub_ogcp1[0] * sub_headgapfreq2 + sub_fgcp1[i - 1] * sub_gapfreq2f[0]);
    }
    if (headgp || jst != 0) {
        for (j = 1; j < sub_lgth2 + 1; j++)
            currentw[j] += (sub_ogcp2[0] * sub_headgapfreq1 + sub_fgcp2[j - 1] * sub_gapfreq1f[0]);
    }
    for (j = 1; j < sub_lgth2 + 1; ++j) {
        m[j] = currentw[j - 1] + sub_ogcp1[1] * sub_gapfreq2f[j - 1];
        mp[j] = 0;
    }
    lastverticalw[0] = currentw[sub_lgth2 - 1];

    int lasti = (tailgp || jen != lgth2 - 1) ? sub_lgth1 + 1 : sub_lgth1;
    int lastj = sub_lgth2 + 1;
    for (i = 1; i < lasti; i++) {
        double *wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;
        previousw[0] = initverticalw[i - 1];
        rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist + i, jst, sub_lgth2, doublework, intwork, 0);
        currentw[0] = initverticalw[i];

        double mi = previousw[0] + sub_ogcp2[1] * sub_gapfreq1f[i - 1];
        int mpi = 0;
        int *ijppt = ijp[i] + 1;
        double *mjpt = m + 1;
        double *prept = previousw;
        double *curpt = currentw + 1;
        int *mpjpt = mp + 1;
        for (j = 1; j < lastj; j++) {
            double wm = *prept;
            *ijppt = 0;
            double g = mi + sub_fgcp2[j - 1] * sub_gapfreq1f[i];
            if (g > wm) { wm = g; *ijppt = -(j - mpi); }
            g = *prept + sub_ogcp2[j] * sub_gapfreq1f[i - 1];
            if (g >= mi) { mi = g; mpi = j - 1; }
            g = *mjpt + sub_fgcp1[i - 1] * sub_gapfreq2f[j];
            if (g > wm) { wm = g; *ijppt = +(i - *mpjpt); }
            g = *prept + sub_ogcp1[i] * sub_gapfreq2f[j - 1];
            if (g >= *mjpt) { *mjpt = g; *mpjpt = i - 1; }
            *curpt += wm;
            ijppt++; mjpt++; prept++; mpjpt++; curpt++;
        }
        lastverticalw[i] = currentw[sub_lgth2 - 1];
    }

    /* Traceback (Atracking). For tailgp=1, just trace from corner. */
    char *gt1 = AllocateCharVec(sub_lgth1 + sub_lgth2 + 3);
    char *gt2 = AllocateCharVec(sub_lgth1 + sub_lgth2 + 3);
    char *gaptable1 = gt1 + sub_lgth1 + sub_lgth2;
    char *gaptable2 = gt2 + sub_lgth1 + sub_lgth2;
    *gaptable1 = 0;
    *gaptable2 = 0;

    int iin = sub_lgth1, jin = sub_lgth2;
    /* Initialise ijp boundaries (Atracking lines 568-575) */
    for (i = 0; i < sub_lgth1 + 1; i++) ijp[i][0] = i + 1;
    for (j = 0; j < sub_lgth2 + 1; j++) ijp[0][j] = -(j + 1);

    /* For non-tailgp: scan last row/col for best terminal */
    if (!tailgp && (ien == lgth1 - 1 || jen == lgth2 - 1)) {
        double wm = lastverticalw[0];
        for (i = 0; i < sub_lgth1; i++) {
            if (lastverticalw[i] >= wm) {
                wm = lastverticalw[i];
                iin = i; jin = sub_lgth2 - 1;
                ijp[sub_lgth1][sub_lgth2] = +(sub_lgth1 - i);
            }
        }
        /* lasthorizontalw = the last row of h, which we don't store
           fully — approximated by currentw (last filled row). */
        for (j = 0; j < sub_lgth2; j++) {
            if (currentw[j] >= wm) {
                wm = currentw[j];
                iin = sub_lgth1 - 1; jin = j;
                ijp[sub_lgth1][sub_lgth2] = -(sub_lgth2 - j);
            }
        }
    }

    int k, klim = sub_lgth1 + sub_lgth2;
    iin = sub_lgth1; jin = sub_lgth2;
    int ifi, jfi, l;
    for (k = 0; k <= klim; k++) {
        if (ijp[iin][jin] < 0) {
            ifi = iin - 1;
            jfi = jin + ijp[iin][jin];
        } else if (ijp[iin][jin] > 0) {
            ifi = iin - ijp[iin][jin];
            jfi = jin - 1;
        } else {
            ifi = iin - 1;
            jfi = jin - 1;
        }
        l = iin - ifi;
        while (--l) {
            *--gaptable1 = 'o';
            *--gaptable2 = '-';
            k++;
        }
        l = jin - jfi;
        while (--l) {
            *--gaptable1 = '-';
            *--gaptable2 = 'o';
            k++;
        }
        if (iin <= 0 || jin <= 0) break;
        *--gaptable1 = 'o';
        *--gaptable2 = 'o';
        k++;
        iin = ifi;
        jin = jfi;
    }

    /* Build out_seq1/out_seq2 from gaptable. */
    int width = 0;
    int p = 0, q = 0;
    char *gptr1 = gaptable1;
    char *gptr2 = gaptable2;
    while (*gptr1 && *gptr2) {
        if (*gptr1 == 'o') {
            out_seq1[width] = seq1[ist + p];
            p++;
        } else {
            out_seq1[width] = '-';
        }
        if (*gptr2 == 'o') {
            out_seq2[width] = seq2[jst + q];
            q++;
        } else {
            out_seq2[width] = '-';
        }
        width++;
        gptr1++;
        gptr2++;
    }
    out_seq1[width] = 0;
    out_seq2[width] = 0;
    *out_width = width;

    /* Free. */
    free(gt1);
    free(gt2);
    FreeFloatMtx(cpmx1);
    FreeFloatMtx(cpmx2);
    FreeFloatVec(ogcp1opt);
    FreeFloatVec(ogcp2opt);
    FreeFloatVec(fgcp1opt);
    FreeFloatVec(fgcp2opt);
    FreeFloatVec(gapfreq1f);
    FreeFloatVec(gapfreq2f);
    FreeFloatVec(ogcp1);
    FreeFloatVec(ogcp2);
    FreeFloatVec(fgcp1);
    FreeFloatVec(fgcp2);
    FreeFloatVec(w1);
    FreeFloatVec(w2);
    FreeFloatVec(initverticalw);
    FreeFloatVec(lastverticalw);
    FreeFloatVec(m);
    FreeIntVec(mp);
    FreeFloatMtx(doublework);
    FreeIntMtx(intwork);
    FreeIntMtx(intmtx);
}

/* =====================================================================
 * Recursive-trace harness: faithful port of MSalignmm_rec that prints
 * level-by-level state to stderr. Lets us see exactly how C composes
 * top + inter-half + bottom into the final alignment for a given
 * input, which is what we need to find the last column of divergence
 * between our Rust `msalignmm` (152 wide on the asymmetric input) and
 * C MAFFT's real `MSalignmm` (151 wide on the same).
 *
 * The trace is gated on `getenv("MSALIGN_TRACE")` so callers don't pay
 * for output unless they want it.
 *
 * To keep this manageable, only single-sequence groups (icyc = jcyc =
 * 1), MEMSAVE=1, no constraints, no warp DP. Same scope as the rest of
 * this file.
 * ===================================================================== */

#define DPTANNI 100

/* Forward decl */
static double rs_msi_rec_traced(
    double **n_dynamicmtx,
    char *seq1, char *seq2,
    double **cpmx1, double **cpmx2,
    int ist, int ien, int jst, int jen,
    int alloclen, int fulllen1, int fulllen2,
    char *mseq1, char *mseq2,
    int depth,
    double **gapinfo,
    int headgp, int tailgp,
    double headgapfreq1_g, double headgapfreq2_g,
    int trace);

/* Base case: faithful copy of MSalignmm_tanni for our scope. Reads
 * gapinfo[0..5] for the gap arrays, slices per (ist, jst). Writes the
 * aligned strings into mseq1 / mseq2 (caller allocates). Returns the
 * alignment score and writes the trace width via the buffer.
 */
static double rs_msi_tanni_inner(
    double **n_dynamicmtx,
    char *seq1, char *seq2,
    double **cpmx1, double **cpmx2,
    int ist, int ien, int jst, int jen,
    int fulllen1, int fulllen2,
    char *mseq1, char *mseq2,
    double **gapinfo,
    int headgp, int tailgp,
    double headgapfreq1_g, double headgapfreq2_g)
{
    int i, j;
    int lgth1 = ien - ist + 1;
    int lgth2 = jen - jst + 1;
    int ll1 = lgth1 + 100, ll2 = lgth2 + 100;
    double wm = 0.0, g;

    double *ogcp1 = gapinfo[0] + ist;
    double *fgcp1 = gapinfo[1] + ist;
    double *ogcp2 = gapinfo[2] + jst;
    double *fgcp2 = gapinfo[3] + jst;
    double *gapfreq1f = gapinfo[4] + ist;
    double *gapfreq2f = gapinfo[5] + jst;
    double headgapfreq1 = (ist > 0) ? gapfreq1f[-1] : headgapfreq1_g;
    double headgapfreq2 = (jst > 0) ? gapfreq2f[-1] : headgapfreq2_g;

    double *w1 = AllocateFloatVec(ll2 + 2);
    double *w2 = AllocateFloatVec(ll2 + 2);
    double *initverticalw = AllocateFloatVec(ll1 + 2);
    double *lastverticalw = AllocateFloatVec(ll1 + 2);
    double *m = AllocateFloatVec(ll2 + 2);
    int *mp = AllocateIntVec(ll2 + 2);
    double **doublework = AllocateFloatMtx(ll1 > ll2 ? ll1 + 2 : ll2 + 2, nalphabets + 1);
    int **intwork = AllocateIntMtx(ll1 > ll2 ? ll1 + 2 : ll2 + 2, nalphabets + 1);
    int **ijp = AllocateIntMtx(ll1 + 1, ll2 + 1);

    double *currentw = w1;
    double *previousw = w2;

    rs_msi_match_calc(n_dynamicmtx, initverticalw, cpmx2, cpmx1, jst, ist, lgth1, doublework, intwork, 1);
    rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist, jst, lgth2, doublework, intwork, 1);

    if (headgp || ist != 0) {
        for (i = 1; i < lgth1 + 1; i++)
            initverticalw[i] += (ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2f[0]);
    }
    if (headgp || jst != 0) {
        for (j = 1; j < lgth2 + 1; j++)
            currentw[j] += (ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1f[0]);
    }
    for (j = 1; j < lgth2 + 1; ++j) {
        m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2f[j - 1];
        mp[j] = 0;
    }
    lastverticalw[0] = currentw[lgth2 - 1];

    int lasti = (tailgp || jen != fulllen2 - 1) ? lgth1 + 1 : lgth1;
    int lastj = lgth2 + 1;
    for (i = 1; i < lasti; i++) {
        double *wtmp = previousw;
        previousw = currentw;
        currentw = wtmp;
        previousw[0] = initverticalw[i - 1];
        rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist + i, jst, lgth2, doublework, intwork, 0);
        currentw[0] = initverticalw[i];

        double mi = previousw[0] + ogcp2[1] * gapfreq1f[i - 1];
        int mpi = 0;
        int *ijppt = ijp[i] + 1;
        double *mjpt = m + 1;
        double *prept = previousw;
        double *curpt = currentw + 1;
        int *mpjpt = mp + 1;
        for (j = 1; j < lastj; j++) {
            wm = *prept;
            *ijppt = 0;
            g = mi + fgcp2[j - 1] * gapfreq1f[i];
            if (g > wm) { wm = g; *ijppt = -(j - mpi); }
            g = *prept + ogcp2[j] * gapfreq1f[i - 1];
            if (g >= mi) { mi = g; mpi = j - 1; }
            g = *mjpt + fgcp1[i - 1] * gapfreq2f[j];
            if (g > wm) { wm = g; *ijppt = +(i - *mpjpt); }
            g = *prept + ogcp1[i] * gapfreq2f[j - 1];
            if (g >= *mjpt) { *mjpt = g; *mpjpt = i - 1; }
            *curpt += wm;
            ijppt++; mjpt++; prept++; mpjpt++; curpt++;
        }
        lastverticalw[i] = currentw[lgth2 - 1];
    }

    /* Atracking */
    char *gt1 = AllocateCharVec(lgth1 + lgth2 + 3);
    char *gt2 = AllocateCharVec(lgth1 + lgth2 + 3);
    char *gaptable1 = gt1 + lgth1 + lgth2;
    char *gaptable2 = gt2 + lgth1 + lgth2;
    *gaptable1 = 0;
    *gaptable2 = 0;

    int iin, jin, ifi, jfi, l, k, klim;
    double *lasthorizontalw = currentw;

    for (i = 0; i < lgth1 + 1; i++) ijp[i][0] = i + 1;
    for (j = 0; j < lgth2 + 1; j++) ijp[0][j] = -(j + 1);

    if (tailgp == 1) {
        /* no special tail handling */
    } else if (ien == fulllen1 - 1 || jen == fulllen2 - 1) {
        wm = lastverticalw[0];
        for (i = 0; i < lgth1; i++) {
            if (lastverticalw[i] >= wm) {
                wm = lastverticalw[i];
                ijp[lgth1][lgth2] = +(lgth1 - i);
            }
        }
        for (j = 0; j < lgth2; j++) {
            if (lasthorizontalw[j] >= wm) {
                wm = lasthorizontalw[j];
                ijp[lgth1][lgth2] = -(lgth2 - j);
            }
        }
    }

    iin = lgth1; jin = lgth2;
    klim = lgth1 + lgth2;
    for (k = 0; k <= klim; k++) {
        if (ijp[iin][jin] < 0) {
            ifi = iin - 1;
            jfi = jin + ijp[iin][jin];
        } else if (ijp[iin][jin] > 0) {
            ifi = iin - ijp[iin][jin];
            jfi = jin - 1;
        } else {
            ifi = iin - 1; jfi = jin - 1;
        }
        l = iin - ifi;
        while (--l) { *--gaptable1 = 'o'; *--gaptable2 = '-'; k++; }
        l = jin - jfi;
        while (--l) { *--gaptable1 = '-'; *--gaptable2 = 'o'; k++; }
        if (iin <= 0 || jin <= 0) break;
        *--gaptable1 = 'o'; *--gaptable2 = 'o'; k++;
        iin = ifi; jin = jfi;
    }

    int width = 0;
    int p = 0, q = 0;
    char *g1 = gaptable1, *g2 = gaptable2;
    while (*g1 && *g2) {
        if (*g1 == 'o') { mseq1[width] = seq1[ist + p]; p++; } else mseq1[width] = '-';
        if (*g2 == 'o') { mseq2[width] = seq2[jst + q]; q++; } else mseq2[width] = '-';
        width++; g1++; g2++;
    }
    mseq1[width] = 0; mseq2[width] = 0;

    free(gt1); free(gt2);
    FreeFloatVec(w1); FreeFloatVec(w2);
    FreeFloatVec(initverticalw); FreeFloatVec(lastverticalw);
    FreeFloatVec(m); FreeIntVec(mp);
    FreeFloatMtx(doublework); FreeIntMtx(intwork);
    FreeIntMtx(ijp);
    return wm;
}

/* Recursive function with trace prints. */
static double rs_msi_rec_traced(
    double **n_dynamicmtx,
    char *seq1, char *seq2,
    double **cpmx1, double **cpmx2,
    int ist, int ien, int jst, int jen,
    int alloclen, int fulllen1, int fulllen2,
    char *mseq1, char *mseq2,
    int depth,
    double **gapinfo,
    int headgp, int tailgp,
    double headgapfreq1_g, double headgapfreq2_g,
    int trace)
{
    int i, j;
    int lgth1 = ien - ist + 1;
    int lgth2 = jen - jst + 1;
    double value = 0.0;
    double wm = 0.0, g;
    int imid = lgth1 / 2;
    int alnlen;

    if (trace) {
        fprintf(stderr, "[TRACE depth=%d] ENTER ist=%d ien=%d jst=%d jen=%d lgth1=%d lgth2=%d headgp=%d tailgp=%d\n",
                depth, ist, ien, jst, jen, lgth1, lgth2, headgp, tailgp);
    }

    if (lgth2 <= 0) {
        for (j = 0; j < lgth1; j++) { mseq1[j] = seq1[ist + j]; mseq2[j] = '-'; }
        mseq1[lgth1] = 0; mseq2[lgth1] = 0;
        if (trace) fprintf(stderr, "[TRACE depth=%d] DEGENERATE return width=%d\n", depth, lgth1);
        return 0.0;
    }

    if (lgth1 < DPTANNI || lgth2 < DPTANNI) {
        value = rs_msi_tanni_inner(n_dynamicmtx, seq1, seq2, cpmx1, cpmx2,
            ist, ien, jst, jen, fulllen1, fulllen2, mseq1, mseq2, gapinfo,
            headgp, tailgp, headgapfreq1_g, headgapfreq2_g);
        if (trace) {
            int w = (int)strlen(mseq1);
            fprintf(stderr, "[TRACE depth=%d] BASE_CASE width=%d score=%.2f\n",
                    depth, w, value);
        }
        return value;
    }

    /* Set up DP scratch (forward + backward + split decision) */
    int ll2 = lgth2 + 100;
    double *w1 = AllocateFloatVec(ll2 + 2);
    double *w2 = AllocateFloatVec(ll2 + 2);
    double *midw = AllocateFloatVec(ll2 + 2);
    double *midm = AllocateFloatVec(ll2 + 2);
    double *midn = AllocateFloatVec(ll2 + 2);
    int *jumpbacki = AllocateIntVec(ll2 + 2);
    int *jumpbackj = AllocateIntVec(ll2 + 2);
    int *jumpforwi = AllocateIntVec(ll2 + 2);
    int *jumpforwj = AllocateIntVec(ll2 + 2);
    int ll1 = lgth1 + 100;
    double *initverticalw = AllocateFloatVec(ll1 + 2);
    double *lastverticalw = AllocateFloatVec(ll1 + 2);
    double *m = AllocateFloatVec(ll2 + 2);
    int *mp = AllocateIntVec(ll2 + 2);
    double **doublework = AllocateFloatMtx(ll1 > ll2 ? ll1 + 2 : ll2 + 2, nalphabets + 1);
    int **intwork = AllocateIntMtx(ll1 > ll2 ? ll1 + 2 : ll2 + 2, nalphabets + 1);

    double *ogcp1 = gapinfo[0] + ist;
    double *fgcp1 = gapinfo[1] + ist;
    double *ogcp2 = gapinfo[2] + jst;
    double *fgcp2 = gapinfo[3] + jst;
    double *gapfreq1f = gapinfo[4] + ist;
    double *gapfreq2f = gapinfo[5] + jst;
    double headgapfreq1 = (ist > 0) ? gapfreq1f[-1] : headgapfreq1_g;
    double headgapfreq2 = (jst > 0) ? gapfreq2f[-1] : headgapfreq2_g;

    double *currentw = w1;
    double *previousw = w2;

    /* Forward DP up to imid */
    rs_msi_match_calc(n_dynamicmtx, initverticalw, cpmx2, cpmx1, jst, ist, lgth1, doublework, intwork, 1);
    rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist, jst, lgth2, doublework, intwork, 1);

    for (i = 1; i < lgth1 + 1; i++)
        initverticalw[i] += (ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2f[0]);
    for (j = 1; j < lgth2 + 1; j++)
        currentw[j] += (ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1f[0]);
    for (j = 1; j < lgth2 + 1; ++j) {
        m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2f[j - 1];
        mp[j] = 0;
    }
    lastverticalw[0] = currentw[lgth2 - 1];

    int jumpi = 0, jumpj = 0, jmid = 0;
    double maxwm = 0.0;

    for (i = 1; i <= imid; i++) {
        double *wtmp = previousw; previousw = currentw; currentw = wtmp;
        previousw[0] = initverticalw[i - 1];
        rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist + i, jst, lgth2, doublework, intwork, 0);
        currentw[0] = initverticalw[i];
        m[0] = ogcp1[i];
        if (i == imid) midm[0] = m[0];

        double mi = previousw[0] + ogcp2[1] * gapfreq1f[i - 1];
        int mpi = 0;
        double *mjpt = m + 1;
        double *prept = previousw;
        double *curpt = currentw + 1;
        int *mpjpt = mp + 1;

        for (j = 1; j < lgth2 + 1; j++) {
            wm = *prept;
            g = mi + fgcp2[j - 1] * gapfreq1f[i];
            if (g > wm) wm = g;
            g = *prept + ogcp2[j] * gapfreq1f[i - 1];
            if (g >= mi) { mi = g; mpi = j - 1; }
            g = *mjpt + fgcp1[i - 1] * gapfreq2f[j];
            if (g > wm) wm = g;
            g = *prept + ogcp1[i] * gapfreq2f[j - 1];
            if (g >= *mjpt) { *mjpt = g; *mpjpt = i - 1; }
            *curpt += wm;

            if (i == imid) {
                jumpbackj[j] = *mpjpt;
                jumpbacki[j] = mpi;
                midw[j] = *curpt;
                midm[j] = *mjpt;
                midn[j] = mi;
            }
            mjpt++; prept++; mpjpt++; curpt++;
        }
        lastverticalw[i] = currentw[lgth2 - 1];
    }

    if (trace) fprintf(stderr, "[MINE] FORWARD_END midw[94]=%.2f midw[95]=%.2f midw[96]=%.2f\n", midw[94], midw[95], midw[96]);

    /* Backward DP */
    rs_msi_match_calc(n_dynamicmtx, initverticalw, cpmx2, cpmx1, jst + lgth2 - 1, ist, lgth1, doublework, intwork, 1);
    rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist + lgth1 - 1, jst, lgth2, doublework, intwork, 1);

    for (i = 0; i < lgth1 - 1; i++)
        initverticalw[i] += (fgcp1[lgth1 - 1] * gapfreq2f[lgth2] + ogcp1[i + 1] * gapfreq2f[lgth2 - 1]);
    for (j = 0; j < lgth2 - 1; j++)
        currentw[j] += (fgcp2[lgth2 - 1] * gapfreq1f[lgth1] + ogcp2[j + 1] * gapfreq1f[lgth1 - 1]);

    for (j = lgth2 - 1; j > -1; --j) {
        m[j] = currentw[j + 1] + fgcp1[lgth1 - 2] * gapfreq2f[j + 1];
        mp[j] = lgth1 - 1;
    }

    double firstm = -9999999.9;
    int firstmp = lgth1;

    for (i = lgth1 - 2; i > -1; i--) {
        double *wtmp = previousw; previousw = currentw; currentw = wtmp;
        previousw[lgth2 - 1] = initverticalw[i + 1];
        rs_msi_match_calc(n_dynamicmtx, currentw, cpmx1, cpmx2, ist + i, jst, lgth2, doublework, intwork, 0);
        currentw[lgth2 - 1] = initverticalw[i];

        double mi = previousw[lgth2 - 1] + fgcp2[lgth2 - 2] * gapfreq1f[i + 1];
        int mpi = lgth2 - 1;
        double *mjpt = m + lgth2 - 2;
        double *prept = previousw + lgth2 - 1;
        double *curpt = currentw + lgth2 - 2;
        int *mpjpt = mp + lgth2 - 2;

        for (j = lgth2 - 2; j > -1; j--) {
            wm = *prept;
            int ijpi = i + 1, ijpj = j + 1;
            g = mi + ogcp2[j + 1] * gapfreq1f[i];
            if (g > wm) { wm = g; ijpj = mpi; ijpi = i + 1; }
            g = *prept + fgcp2[j] * gapfreq1f[i + 1];
            if (g >= mi) { mi = g; mpi = j + 1; }
            g = *mjpt + ogcp1[i + 1] * gapfreq2f[j];
            if (g > wm) { wm = g; ijpi = *mpjpt; ijpj = j + 1; }
            g = *prept + fgcp1[i] * gapfreq2f[j + 1];
            if (g >= *mjpt) { *mjpt = g; *mpjpt = i + 1; }

            if (i == jumpi || i == imid - 1) {
                jumpforwi[j] = ijpi;
                jumpforwj[j] = ijpj;
            }
            if (i == imid) { midw[j] += wm; midm[j + 1] += *mjpt; }
            if (i == imid - 1) midn[j] += mi;
            *curpt += wm;
            mjpt--; prept--; mpjpt--; curpt--;
        }
        g = *prept + fgcp1[i];
        if (firstm < g) { firstm = g; firstmp = i + 1; }
        if (i == imid) midm[j + 1] += firstm;

        if (i == imid - 1) {
            if (trace) fprintf(stderr, "[MINE] midw[94]=%.2f midw[95]=%.2f midw[96]=%.2f midw[99]=%.2f\n", midw[94], midw[95], midw[96], midw[99]);
            if (trace) fprintf(stderr, "[MINE] midm[94]=%.2f midm[95]=%.2f midm[96]=%.2f\n", midm[94], midm[95], midm[96]);
            if (trace) fprintf(stderr, "[MINE] midn[93]=%.2f midn[94]=%.2f midn[95]=%.2f\n", midn[93], midn[94], midn[95]);
            maxwm = midw[1]; jmid = 0;
            for (j = 2; j < lgth2 - 1; j++) {
                if (midw[j] > maxwm) { jmid = j; maxwm = midw[j]; }
            }
            for (j = 0; j < lgth2 + 1; j++) {
                if (midm[j] > maxwm) { jmid = j; maxwm = midm[j]; }
            }
            wm = midw[jmid]; jumpi = imid - 1; jumpj = jmid - 1;
            if (jmid > 0 && midn[jmid - 1] > wm) {
                jumpj = jumpbacki[jmid]; wm = midn[jmid - 1];
            }
            if (midm[jmid] > wm) {
                jumpi = jumpbackj[jmid]; jumpj = jmid - 1; wm = midm[jmid];
            }
        }
        if (i == jumpi) {
            if (jmid == 0) {
                jumpj = 0; jmid = 1;
                if (imid < firstmp - 1) { jumpi = firstmp; imid = firstmp + 1; }
            } else if (jmid >= lgth2) {
                jumpi = imid - 1; jmid = lgth2; jumpj = lgth2 - 1;
            } else {
                imid = jumpforwi[jumpj];
                jmid = jumpforwj[jumpj];
                if (imid == jumpi) jumpi = imid - 1;
            }
            break;
        }
    }

    if (trace) {
        fprintf(stderr, "[TRACE depth=%d] SPLIT imid=%d jmid=%d jumpi=%d jumpj=%d\n",
                depth, imid, jmid, jumpi, jumpj);
    }

    FreeFloatVec(w1); FreeFloatVec(w2);
    FreeFloatVec(midw); FreeFloatVec(midm); FreeFloatVec(midn);
    FreeIntVec(jumpbacki); FreeIntVec(jumpbackj);
    FreeIntVec(jumpforwi); FreeIntVec(jumpforwj);
    FreeFloatVec(initverticalw); FreeFloatVec(lastverticalw);
    FreeFloatVec(m); FreeIntVec(mp);
    FreeFloatMtx(doublework); FreeIntMtx(intwork);

    /* Top recursion */
    char *aseq1 = mseq1;
    char *aseq2 = mseq2;
    value = rs_msi_rec_traced(n_dynamicmtx, seq1, seq2, cpmx1, cpmx2,
        ist, ist + jumpi, jst, jst + jumpj, alloclen,
        fulllen1, fulllen2, aseq1, aseq2, depth + 1, gapinfo,
        headgp, tailgp, headgapfreq1_g, headgapfreq2_g, trace);
    if (trace) {
        fprintf(stderr, "[TRACE depth=%d] TOP_DONE strlen(mseq1)=%d\n",
                depth, (int)strlen(mseq1));
    }

    /* Inter-half horizontal */
    int len = (int)strlen(mseq1);
    int l_horiz = jmid - jumpj - 1;
    if (l_horiz > 0) {
        for (i = 0; i < l_horiz; i++) {
            mseq1[len + i] = '-';
            mseq2[len + i] = seq2[jst + jumpj + 1 + i];
        }
        mseq1[len + l_horiz] = 0;
        mseq2[len + l_horiz] = 0;
        value += (ogcp2[jumpj + 1] + fgcp2[jmid - 1]);
        if (trace) {
            fprintf(stderr, "[TRACE depth=%d] INTER_HORIZ l=%d strlen=%d\n",
                    depth, l_horiz, (int)strlen(mseq1));
        }
    }

    /* Inter-half vertical */
    len = (int)strlen(mseq1);
    int l_vert = imid - jumpi - 1;
    if (l_vert > 0) {
        for (i = 0; i < l_vert; i++) {
            mseq1[len + i] = seq1[ist + jumpi + 1 + i];
            mseq2[len + i] = '-';
        }
        mseq1[len + l_vert] = 0;
        mseq2[len + l_vert] = 0;
        value += (ogcp1[jumpi + 1] + fgcp1[imid - 1]);
        if (trace) {
            fprintf(stderr, "[TRACE depth=%d] INTER_VERT l=%d strlen=%d\n",
                    depth, l_vert, (int)strlen(mseq1));
        }
    }

    /* Advance pointers (MEMSAVE) */
    alnlen = (int)strlen(mseq1);
    aseq1 += alnlen;
    aseq2 += alnlen;

    /* Bottom recursion */
    if (trace) {
        fprintf(stderr, "[TRACE depth=%d] BOTTOM_RECURSE_PRE ist=%d ien=%d jst=%d jen=%d (offset=%d)\n",
                depth, ist + imid, ien, jst + jmid, jen, alnlen);
    }
    value += rs_msi_rec_traced(n_dynamicmtx, seq1, seq2, cpmx1, cpmx2,
        ist + imid, ien, jst + jmid, jen, alloclen,
        fulllen1, fulllen2, aseq1, aseq2, depth + 1, gapinfo,
        headgp, tailgp, headgapfreq1_g, headgapfreq2_g, trace);
    if (trace) {
        fprintf(stderr, "[TRACE depth=%d] BOTTOM_DONE total_strlen=%d\n",
                depth, (int)strlen(mseq1));
    }

    return value;
}

/*
 * Public entry point. Sets up cpmx + gapinfo for the FULL inputs,
 * then runs `rs_msi_rec_traced` (faithful C re-implementation of
 * MSalignmm_rec) with the trace flag. Returns the aligned strings
 * (caller allocates) and width.
 *
 * If `getenv("MSALIGN_TRACE")` is set the recursion prints
 * level-by-level state to stderr.
 */
void rs_msalignmm_full_trace(
    double **n_dynamicmtx,
    char *seq1, char *seq2,
    int lgth1, int lgth2,
    int headgp, int tailgp,
    char *out_seq1, char *out_seq2, int *out_width)
{
    int i;
    int ll1 = lgth1 + 100, ll2 = lgth2 + 100;
    char *seqs1[1] = {seq1};
    char *seqs2[1] = {seq2};
    double eff1[1] = {1.0};
    double eff2[1] = {1.0};

    double **cpmx1 = AllocateFloatMtx(nalphabets, ll1 + 2);
    double **cpmx2 = AllocateFloatMtx(nalphabets, ll2 + 2);
    cpmx_calc_new(seqs1, cpmx1, eff1, lgth1, 1);
    cpmx_calc_new(seqs2, cpmx2, eff2, lgth2, 1);

    double *ogcp1opt = AllocateFloatVec(ll1 + 2);
    double *ogcp2opt = AllocateFloatVec(ll2 + 2);
    double *fgcp1opt = AllocateFloatVec(ll1 + 2);
    double *fgcp2opt = AllocateFloatVec(ll2 + 2);
    st_OpeningGapCount(ogcp1opt, 1, seqs1, eff1, lgth1);
    st_FinalGapCount(fgcp1opt, 1, seqs1, eff1, lgth1);
    st_OpeningGapCount(ogcp2opt, 1, seqs2, eff2, lgth2);
    st_FinalGapCount(fgcp2opt, 1, seqs2, eff2, lgth2);

    double *gapfreq1f = AllocateFloatVec(ll1 + 2);
    double *gapfreq2f = AllocateFloatVec(ll2 + 2);
    gapcountf(gapfreq1f, seqs1, 1, eff1, lgth1);
    gapcountf(gapfreq2f, seqs2, 1, eff2, lgth2);
    for (i = 0; i < lgth1 + 1; i++) gapfreq1f[i] = 1.0 - gapfreq1f[i];
    for (i = 0; i < lgth2 + 1; i++) gapfreq2f[i] = 1.0 - gapfreq2f[i];
    double headgapfreq1 = 1.0, headgapfreq2 = 1.0;

    double fpenalty = (double)penalty;
    double *ogcp1 = AllocateFloatVec(ll1 + 2);
    double *ogcp2 = AllocateFloatVec(ll2 + 2);
    double *fgcp1 = AllocateFloatVec(ll1 + 2);
    double *fgcp2 = AllocateFloatVec(ll2 + 2);
    for (i = 0; i < lgth1; i++) {
        ogcp1[i] = 0.5 * (1.0 - ogcp1opt[i]) * fpenalty * gapfreq1f[i];
        fgcp1[i] = 0.5 * (1.0 - fgcp1opt[i]) * fpenalty * gapfreq1f[i];
    }
    for (i = 0; i < lgth2; i++) {
        ogcp2[i] = 0.5 * (1.0 - ogcp2opt[i]) * fpenalty * gapfreq2f[i];
        fgcp2[i] = 0.5 * (1.0 - fgcp2opt[i]) * fpenalty * gapfreq2f[i];
    }

    double *gapinfo[6] = {ogcp1, fgcp1, ogcp2, fgcp2, gapfreq1f, gapfreq2f};

    int trace = (getenv("MSALIGN_TRACE") != NULL) ? 1 : 0;
    out_seq1[0] = 0;
    out_seq2[0] = 0;

    rs_msi_rec_traced(n_dynamicmtx, seq1, seq2, cpmx1, cpmx2,
        0, lgth1 - 1, 0, lgth2 - 1, lgth1 + lgth2 + 100,
        lgth1, lgth2, out_seq1, out_seq2, 0, gapinfo,
        headgp, tailgp, headgapfreq1, headgapfreq2, trace);

    *out_width = (int)strlen(out_seq1);

    FreeFloatMtx(cpmx1); FreeFloatMtx(cpmx2);
    FreeFloatVec(ogcp1opt); FreeFloatVec(ogcp2opt);
    FreeFloatVec(fgcp1opt); FreeFloatVec(fgcp2opt);
    FreeFloatVec(gapfreq1f); FreeFloatVec(gapfreq2f);
    FreeFloatVec(ogcp1); FreeFloatVec(ogcp2);
    FreeFloatVec(fgcp1); FreeFloatVec(fgcp2);
}
