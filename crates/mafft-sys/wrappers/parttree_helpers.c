/// Re-exports for `splittbfast`/`disttbfast` 6-mer helpers that are
/// otherwise stuck behind `static` (splittbfast.c) or live in C files with
/// their own `main()` (disttbfast.c, addsingle.c). The bodies are copies of
/// the canonical versions in mltaln9.c / addsingle.c so we can FFI-validate
/// the Rust port without dragging a `main` symbol into the test binary.

#include "mltaln.h"
#include <stdio.h>
#include <stdlib.h>

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
