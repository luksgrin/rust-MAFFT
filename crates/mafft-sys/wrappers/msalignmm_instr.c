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
