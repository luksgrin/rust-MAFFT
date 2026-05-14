/// Re-exports for `splittbfast`/`disttbfast` 6-mer helpers that are
/// otherwise stuck behind `static` (splittbfast.c) or live in C files with
/// their own `main()` (disttbfast.c, addsingle.c). The bodies are copies of
/// the canonical versions in mltaln9.c / addsingle.c so we can FFI-validate
/// the Rust port without dragging a `main` symbol into the test binary.

#include "mltaln.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char amino_grp[0x100];

/// `nunknown` is `static` in the upstream C files (addsingle.c:806,
/// disttbfast.c:683) so we provide a local stand-in for the wrapper.
static int nunknown;

/// Copy of `addsingle.c::seq_grp` — protein 6-group encoding,
/// drops chars whose `amino_grp[c] >= 6` (X, '.', '-').
int seq_grp( int *grp, char *seq )
{
    int tmp;
    int *grpbk = grp;
    while( *seq )
    {
        tmp = amino_grp[(int)(unsigned char)*seq++];
        if( tmp < 6 )
            *grp++ = tmp;
        else
            nunknown++;
    }
    *grp = -1; // END_OF_VEC
    return( grp - grpbk );
}

/// Copy of `addsingle.c::seq_grp_nuc` — DNA 4-group encoding.
int seq_grp_nuc( int *grp, char *seq )
{
    int tmp;
    int *grpbk = grp;
    while( *seq )
    {
        tmp = amino_grp[(int)(unsigned char)*seq++];
        if( tmp < 4 )
            *grp++ = tmp;
        else
            nunknown++;
    }
    *grp = -1;
    return( grp - grpbk );
}

/// Copy of `addsingle.c::makecompositiontable_p` — 6-mer frequency table.
void makecompositiontable_p( int *table, int *pointt )
{
    int point;
    while( ( point = *pointt++ ) != -1 )
    {
        table[point]++;
    }
}

/// Copy of `addsingle.c::makepointtable` — protein rolling 6-mer encoder.
void makepointtable( int *pointt, int *n )
{
    int point;
    int *p;

    p = n;
    point  = *n++ *  7776;
    point += *n++ *  1296;
    point += *n++ *   216;
    point += *n++ *    36;
    point += *n++ *     6;
    point += *n++;
    *pointt++ = point;

    while( *n != -1 )
    {
        point -= *p++ * 7776;
        point *= 6;
        point += *n++;
        *pointt++ = point;
    }
    *pointt = -1;
}

/// Copy of `addsingle.c::makepointtable_nuc` — DNA rolling 6-mer encoder.
void makepointtable_nuc( int *pointt, int *n )
{
    int point;
    int *p;

    p = n;
    point  = *n++ *  1024;
    point += *n++ *   256;
    point += *n++ *    64;
    point += *n++ *    16;
    point += *n++ *     4;
    point += *n++;
    *pointt++ = point;

    while( *n != -1 )
    {
        point -= *p++ * 1024;
        point *= 4;
        point += *n++;
        *pointt++ = point;
    }
    *pointt = -1;
}

/// Copy of `disttbfast.c::preferenceval` — 1e-14 tie-breaker for the
/// initial pairwise scan in compactdisthalfmtxthread.
static double rs_preferenceval( int ori, int pos, int max )
{
    pos -= ori;
    if( pos < 0 ) pos += max;
    return( 0.00000000000001 * pos );
}

/// Allocate `partmtx[njob]` initialized to NULL (mirrors C
/// `disttbfast.c::preparepartmtx`).
double **rs_prepare_partmtx( int njob )
{
    int i;
    double **m = (double **)calloc( njob, sizeof( double * ) );
    for( i=0; i<njob; i++ ) m[i] = NULL;
    return m;
}

extern void compacttreegivendist(
    int njob, double *mindists, int *neighbors,
    int ***topol, double **len, char **name,
    Treedep *dep, int treeout );

/// Drive C's `compacttreegivendist` (`mltaln9.c:5221`) — the algorithm
/// `--memsavetree` actually uses (NOT `compacttree_memsaveselectable`).
/// Builds a tree from the precomputed `mindist`/`nearest` arrays via
/// stepwise insertion, walking up the tree from each leaf's nearest
/// neighbor.
///
/// Returns the per-step `topol[k][0][0]`/`[1][0]` (rep0/rep1 of the
/// merged subtrees) and `len[k][0]`/`[1]` (branch lengths) via the
/// `out_*` buffers (each length `nseq - 1`).
void rs_compacttreegivendist(
    int nseq,
    double *mindist_in,
    int *nearest_in,
    int *out_topol0,
    int *out_topol1,
    double *out_len0,
    double *out_len1 )
{
    int i;
    int ***topol = (int ***)calloc( nseq, sizeof( int ** ) );
    double **len = (double **)calloc( nseq, sizeof( double * ) );
    Treedep *dep = (Treedep *)calloc( nseq, sizeof( Treedep ) );
    for( i=0; i<nseq; i++ )
    {
        topol[i] = (int **)calloc( 2, sizeof( int * ) );
        len[i] = (double *)calloc( 2, sizeof( double ) );
    }
    double *mindist_work = (double *)malloc( nseq * sizeof( double ) );
    int *nearest_work = (int *)malloc( nseq * sizeof( int ) );
    memcpy( mindist_work, mindist_in, nseq * sizeof( double ) );
    memcpy( nearest_work, nearest_in, nseq * sizeof( int ) );

    njob = nseq;

    compacttreegivendist( nseq, mindist_work, nearest_work, topol, len, NULL, dep, 0 );

    for( i=0; i<nseq-1; i++ )
    {
        out_topol0[i] = topol[i][0] ? topol[i][0][0] : -1;
        out_topol1[i] = topol[i][1] ? topol[i][1][0] : -1;
        out_len0[i] = len[i][0];
        out_len1[i] = len[i][1];
    }

    for( i=0; i<nseq; i++ )
    {
        if( topol[i][0] ) free( topol[i][0] );
        if( topol[i][1] ) free( topol[i][1] );
        free( topol[i] );
        free( len[i] );
    }
    free( topol );
    free( len );
    free( dep );
    free( mindist_work );
    free( nearest_work );
}

/// Drive C's `compacttree_memsaveselectable` with `howcompact=2`,
/// `memsave=1`, `seq=NULL`, `skiptable=NULL`, `dep=NULL`,
/// `treeout=0` — kept for reference but NOT the algorithm
/// `--memsavetree` actually uses (that's `rs_compacttreegivendist`).
void rs_compacttree_memsaveselectable_kmer(
    int nseq,
    int **pointt,
    int *nogaplen,
    int *selfscore,
    double *mindist_in,    // copied in; algorithm mutates internally
    int *nearest_in,
    int *out_topol0,
    int *out_topol1,
    double *out_len0,
    double *out_len1 )
{
    int i;
    double **partmtx = rs_prepare_partmtx( nseq );
    int ***topol = (int ***)calloc( nseq, sizeof( int ** ) );
    double **len = (double **)calloc( nseq, sizeof( double * ) );
    for( i=0; i<nseq; i++ )
    {
        topol[i] = (int **)calloc( 2, sizeof( int * ) );
        // topol[i][0/1] = NULL initially; compacttree_memsaveselectable
        // does its own realloc inside the loop.
        len[i] = (double *)calloc( 2, sizeof( double ) );
    }
    // Working copies of mindist/nearest (algorithm overwrites them).
    double *mindist_work = (double *)malloc( nseq * sizeof( double ) );
    int *nearest_work = (int *)malloc( nseq * sizeof( int ) );
    memcpy( mindist_work, mindist_in, nseq * sizeof( double ) );
    memcpy( nearest_work, nearest_in, nseq * sizeof( int ) );

    // Set globals the algorithm reads beyond its function-parameter
    // shadow:
    //   - `treemethod = 'X'` selects cluster_mix_double linkage.
    //   - `nthreadpair = 0` forces the serial dispatch (the threading
    //     setup needs more globals we don't replicate).
    //   - `njob = nseq` — `compacttree_memsaveselectable` does
    //     `joblist = calloc(njob, sizeof(int))` with the GLOBAL `njob`
    //     (not the function arg). Forgetting this gives joblist size 0
    //     → OOB write in the join-list build → SIGSEGV/SIGBUS.
    treemethod = 'X';
    nthreadpair = 0;
    njob = nseq;

    compacttree_memsaveselectable(
        nseq,
        partmtx,
        nearest_work,
        mindist_work,
        pointt,
        selfscore,
        NULL,        // seq → triggers verycompactkmerdistarrthreadjoblist
        NULL,        // skiptable
        topol,
        len,
        NULL,        // name
        nogaplen,    // nlen
        NULL,        // dep
        0,           // treeout
        2,           // howcompact
        1            // memsave
    );

    // memsave=1 stores [smallest_leaf, -1] in each topol[k][side].
    for( i=0; i<nseq-1; i++ )
    {
        out_topol0[i] = topol[i][0] ? topol[i][0][0] : -1;
        out_topol1[i] = topol[i][1] ? topol[i][1][0] : -1;
        out_len0[i] = len[i][0];
        out_len1[i] = len[i][1];
    }

    // Cleanup.
    for( i=0; i<nseq; i++ )
    {
        if( topol[i][0] ) free( topol[i][0] );
        if( topol[i][1] ) free( topol[i][1] );
        free( topol[i] );
        free( len[i] );
    }
    free( topol );
    free( len );
    free( partmtx );
    free( mindist_work );
    free( nearest_work );

    commonsextet_p( NULL, NULL );
}

/// Re-export of `disttbfast.c::compactdisthalfmtxthread` (static; we
/// copy its loop body here so the Rust FFI test can compare its output
/// against our `initial_mindist` cell-by-cell). Single-threaded only.
///
/// Reads `pointt[]`, `nogaplen[]`, `selfscore[]` (length `nseq`); fills
/// `mindist[nseq]` and `mindistfrom[nseq]` with the smallest-pair
/// distance from each `i` to any `j < i`. The trailing `preferenceval`
/// offset is SUBTRACTED back out at the end (matching
/// `disttbfast.c:3718,3764`).
///
/// `distcompact` and `commonsextet_p` are called from `mltaln9.c` — the
/// caller MUST have initialized C globals (`tsize`, `lenfac{a,b,c,d}`,
/// `maxl`) and called `constants()` before invoking this.
void rs_compact_initial_mindist(
    int nseq,
    int **pointt,
    int *nogaplen,
    int *selfscore,
    double *mindist,
    int *mindistfrom )
{
    int i, j;
    int *table1;

    for( i=0; i<nseq; i++ )
    {
        mindist[i] = 999.9;
        mindistfrom[i] = -1;
    }

    for( i=nseq-1; i>=0; i-- )
    {
        table1 = (int *)calloc( tsize, sizeof( int ) );
        if( !table1 ) { fprintf( stderr, "rs_compact_initial_mindist: calloc\n" ); exit( 1 ); }
        makecompositiontable_p( table1, pointt[i] );

        for( j=i-1; j>-1; j-- )
        {
            double tmpdist = distcompact( nogaplen[i], nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j] );
            double preference = rs_preferenceval( i, j, nseq );
            double tmpdistx = tmpdist + preference;
            if( tmpdistx < mindist[i] )
            {
                mindist[i] = tmpdistx;
                mindistfrom[i] = j;
            }
        }
        free( table1 );
    }
    commonsextet_p( NULL, NULL ); // free static memo

    for( i=0; i<nseq; i++ )
    {
        if( mindistfrom[i] >= 0 )
            mindist[i] -= rs_preferenceval( i, mindistfrom[i], nseq );
    }
}
