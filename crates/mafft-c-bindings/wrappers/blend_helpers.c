/// FFI wrappers for `static` blend helpers in Salignmm.c (createcpmxresult,
/// createogresult, createfgresult, creategapfreqresult). These are called
/// during the cached-profile path in C's progressive alignment but aren't
/// linkable externally; we copy the function bodies verbatim from
/// `mafft-upstream/core/Salignmm.c:608-687,704-763,775-823` so the Rust
/// FFI test can compare our `blend_profiles_exact` cell-by-cell.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int nalphabets;

/// Verbatim copy of `Salignmm.c::createcpmxresult` (lines 608-661).
///
/// Inputs:
/// - `cpmxresult`: array of nalphabets pointers; each filled by this fn.
/// - `limk`: width hint passed to calloc (we pass `alen` from Rust).
/// - `eff1, eff2`: normalized cluster weights.
/// - `cpmx1, cpmx2`: per-cluster cpmx arrays (nalphabets x prof_len).
/// - `gaptable1, gaptable2`: null-terminated `o`/`-` strings of length alen.
///
/// Note: `usehist1/usehist2` are passed 0 in our wrapper — we never want C
/// to free our Rust-allocated children. The free branches in the original
/// are gated behind `usehist1/2` so passing 0 skips them safely.
void rs_createcpmxresult(
    double **cpmxresult,
    int limk,
    double eff1, double eff2,
    double ***cpmx1, double ***cpmx2,
    char *gaptable1, char *gaptable2
)
{
    int i, j, p;
    int alen = strlen( gaptable1 );

    for( i=0; i<nalphabets; i++ )
    {
        cpmxresult[i] = calloc( limk+1, sizeof( double ) );
        for( j=0; j<alen; j++ ) cpmxresult[i][j] = 0.0;
        for( j=0,p=0; j<alen; j++ )
        {
            if( gaptable1[j] == '-' )
                ;
            else
                cpmxresult[i][j] += (*cpmx1)[i][p++] * eff1;
        }

        for( j=0,p=0; j<alen; j++ )
        {
            if( gaptable2[j] == '-' )
                ;
            else
                cpmxresult[i][j] += (*cpmx2)[i][p++] * eff2;
        }
    }
}

/// Verbatim copy of `Salignmm.c::creategapfreqresult` (lines 664-690).
void rs_creategapfreqresult(
    double **gapfresult,
    int limk,
    double eff1, double eff2,
    double *gapf1, double *gapf2,
    char *gaptable1, char *gaptable2
)
{
    int j, p;
    int alen = strlen( gaptable1 );
    (*gapfresult) = calloc( limk+1, sizeof( double ) );

    for( j=0; j<alen+1; j++ ) (*gapfresult)[j] = 0.0;
    for( j=0,p=0; j<alen+1; j++ )
    {
        if( gaptable1[j] == '-' )
            ;
        else
            (*gapfresult)[j] += gapf1[p++] * eff1;
    }

    for( j=0,p=0; j<alen; j++ )
    {
        if( gaptable2[j] == '-' )
            ;
        else
            (*gapfresult)[j] += gapf2[p++] * eff2;
    }
    (*gapfresult)[j] = 1.0;
}

/// Verbatim copy of `Salignmm.c::createogresult` (lines 704-735).
void rs_createogresult(
    double **gapfresult,
    int limk,
    double eff1, double eff2,
    double *ori1, double *ori2,
    double *gf1, double *gf2,
    char *gaptable1, char *gaptable2
)
{
    int j, p;
    int alen = strlen( gaptable1 );

    *gapfresult = calloc( limk+1, sizeof( double ) );

    for( j=0; j<alen; j++ ) (*gapfresult)[j] = 0.0;
    for( j=0,p=0; j<alen; j++ )
    {
        if( gaptable1[j] == '-' )
        {
            if( j==0 )
            {
                (*gapfresult)[j] += 1.0 * eff1;
            }
            else if ( j && gaptable1[j-1] != '-' )
            {
                (*gapfresult)[j] += (gf1[p-1]) * eff1;
            }
        }
        else
        {
            if( j==0 || ( j && gaptable1[j-1] != '-' ) )
            {
                (*gapfresult)[j] += ori1[p] * eff1;
            }
            p++;
        }
    }

    for( j=0,p=0; j<alen; j++ )
    {
        if( gaptable2[j] == '-' )
        {
            if( j==0 )
            {
                (*gapfresult)[j] += 1.0 * eff2;
            }
            else if ( j && gaptable2[j-1] != '-' )
            {
                (*gapfresult)[j] += (gf2[p-1]) * eff2;
            }
        }
        else
        {
            if( j==0 || ( j && gaptable2[j-1] != '-' ) )
            {
                (*gapfresult)[j] += ori2[p] * eff2;
            }
            p++;
        }
    }
}

/// Verbatim copy of `Salignmm.c::createfgresult` (lines 775-825).
void rs_createfgresult(
    double **gapfresult,
    int limk,
    double eff1, double eff2,
    double *ori1, double *ori2,
    double *gf1, double *gf2,
    char *gaptable1, char *gaptable2
)
{
    int j, p;
    int alen = strlen( gaptable1 );
    (*gapfresult) = calloc( limk+1, sizeof( double ) );

    for( j=0; j<alen; j++ ) (*gapfresult)[j] = 0.0;
    for( j=0,p=0; j<alen; j++ )
    {
        if( gaptable1[j] == '-' )
        {
            if( j==alen-1 )
            {
                (*gapfresult)[j] += eff1;
            }
            else if( gaptable1[j+1] != '-' )
            {
                (*gapfresult)[j] += gf1[p] * eff1;
            }
        }
        else
        {
            if( gaptable1[j+1] != '-' )
                (*gapfresult)[j] += ori1[p] * eff1;
            p++;
        }
    }

    for( j=0,p=0; j<alen; j++ )
    {
        if( gaptable2[j] == '-' )
        {
            if( j==alen-1 )
            {
                (*gapfresult)[j] += eff2;
            }
            else if( gaptable2[j+1] != '-' )
            {
                (*gapfresult)[j] += gf2[p] * eff2;
            }
        }
        else
        {
            if( gaptable2[j+1] != '-' )
                (*gapfresult)[j] += ori2[p] * eff2;
            p++;
        }
    }
}
