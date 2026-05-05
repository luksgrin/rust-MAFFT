//! Raw FFI bindings to the MAFFT C library.
//!
//! This crate compiles the MAFFT C sources into a static library and exposes
//! the C function signatures to Rust. It is intended as a transitional crate:
//! as modules are rewritten in Rust, the corresponding C sources will be
//! removed from the build and the FFI declarations replaced with native Rust.

#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
#![allow(clippy::upper_case_acronyms)]

use std::os::raw::{c_char, c_double, c_int, c_void};

// ---------------------------------------------------------------------------
// Core data structures (C layout)
// ---------------------------------------------------------------------------

/// Local homology region between a pair of sequences.
/// Linked list in C; corresponds to `struct _LocalHom` in mltaln.h.
#[repr(C)]
pub struct LocalHom {
    pub next: *mut LocalHom,
    pub last: *mut LocalHom,
    pub start1: c_int,
    pub end1: c_int,
    pub start2: c_int,
    pub end2: c_int,
    pub opt: c_double,
    pub overlapaa: c_int,
    pub extended: c_int,
    pub importance: c_double,
    pub rimportance: c_double,
    pub korh: c_char,
    pub nokori: c_int,
}

/// Guide tree node. Corresponds to `struct _Node` in mltaln.h.
#[repr(C)]
pub struct Node {
    pub children: [*mut Node; 3],
    pub tmpChildren: [c_int; 3],
    pub length: [c_double; 3],
    pub weightptr: [*mut c_double; 3],
    pub top: [c_int; 3],
    pub members: [*mut c_int; 3],
}

/// Aligned segment region. Corresponds to `struct _Segment` in mltaln.h.
#[repr(C)]
pub struct Segment {
    pub start: c_int,
    pub end: c_int,
    pub center: c_int,
    pub score: c_double,
    pub skipForeward: c_int,
    pub skipBackward: c_int,
    pub pair: *mut Segment,
    pub number: c_int,
}

/// Complex number for FFT. Corresponds to `struct _Fukusosuu` in mltaln.h.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Fukusosuu {
    pub R: c_double,
    pub I: c_double,
}

/// Tree dependency info. Corresponds to `struct _Treedep` in mltaln.h.
#[repr(C)]
pub struct Treedep {
    pub child0: c_int,
    pub child1: c_int,
    pub done: c_int,
    pub distfromtip: c_double,
}

/// RNA base pair info. Corresponds to `struct _RNApair` in mltaln.h.
#[repr(C)]
pub struct RNApair {
    pub uppos: c_int,
    pub upscore: c_double,
    pub downpos: c_int,
    pub downscore: c_double,
    pub bestpos: c_int,
    pub bestscore: c_double,
}

/// Addtree info. Corresponds to `struct _Addtree` in mltaln.h.
#[repr(C)]
pub struct Addtree {
    pub nearest: c_int,
    pub dist1: c_double,
    pub neighbors: *mut c_char,
    pub dist2: c_double,
}

/// Gap position. Corresponds to `struct _GapPos` in mltaln.h.
#[repr(C)]
pub struct GapPos {
    pub pos: c_int,
    pub len: c_int,
}

/// Gap pattern. Corresponds to `struct _Gappattern` in mltaln.h.
#[repr(C)]
pub struct Gappat {
    pub len: c_int,
    pub freq: c_double,
}

/// NodeInCub. Corresponds to `struct _NodeInCub` in mltaln.h.
#[repr(C)]
pub struct NodeInCub {
    pub step: c_int,
    pub LorR: c_int,
}

/// Extended anchor. Corresponds to `struct _extanch` in mltaln.h.
#[repr(C)]
pub struct ExtAnch {
    pub i: c_int,
    pub j: c_int,
    pub starti: c_int,
    pub endi: c_int,
    pub startj: c_int,
    pub endj: c_int,
    pub score: c_int,
}

// ---------------------------------------------------------------------------
// Constants from mltaln.h
// ---------------------------------------------------------------------------

pub const M: usize = 500_000;
pub const B: usize = 256;
pub const FFT_THRESHOLD: c_int = 80;
pub const FFT_WINSIZE_P: c_int = 20;
pub const FFT_WINSIZE_D: c_int = 100;
pub const MAXITERATION: c_int = 500;
pub const NOTSPECIFIED: c_int = 100_009;

// ---------------------------------------------------------------------------
// Extern C functions
// ---------------------------------------------------------------------------

unsafe extern "C" {
    // -- Initialization & globals --
    pub fn constants(nseq: c_int, seq: *mut *mut c_char);
    pub fn freeconstants();
    pub fn initglobalvariables();
    pub fn initFiles();
    pub fn closeFiles();

    // -- I/O --
    pub fn PreRead(fp: *mut c_void, locnjob: *mut c_int, locnlenmax: *mut c_int);
    pub fn readData(fp: *mut c_void, name: *mut c_void, nlen: *mut c_int, seq: *mut *mut c_char);
    pub fn readData_pointer(
        fp: *mut c_void,
        name: *mut *mut c_char,
        nlen: *mut c_int,
        seq: *mut *mut c_char,
    );
    pub fn writeData_pointer(
        fp: *mut c_void,
        locnjob: c_int,
        name: *mut *mut c_char,
        nlen: *mut c_int,
        aseq: *mut *mut c_char,
    );

    // -- Matrix allocation (mtxutl) --
    pub fn AllocateCharVec(l1: c_int) -> *mut c_char;
    pub fn AllocateCharMtx(l1: c_int, l2: c_int) -> *mut *mut c_char;
    pub fn FreeCharMtx(mtx: *mut *mut c_char);
    pub fn AllocateFloatVec(l1: c_int) -> *mut c_double;
    pub fn FreeFloatVec(vec: *mut c_double);
    pub fn AllocateFloatMtx(l1: c_int, l2: c_int) -> *mut *mut c_double;
    pub fn FreeFloatMtx(mtx: *mut *mut c_double);
    pub fn AllocateIntVec(l1: c_int) -> *mut c_int;
    pub fn FreeIntVec(vec: *mut c_int);
    pub fn AllocateIntMtx(l1: c_int, l2: c_int) -> *mut *mut c_int;
    pub fn FreeIntMtx(mtx: *mut *mut c_int);
    pub fn AllocateDoubleMtx(l1: c_int, l2: c_int) -> *mut *mut c_double;
    pub fn FreeDoubleMtx(mtx: *mut *mut c_double);
    pub fn AllocateFloatHalfMtx(l1: c_int) -> *mut *mut c_double;
    pub fn FreeFloatHalfMtx(mtx: *mut *mut c_double, n: c_int);

    // -- Scoring & substitution --
    pub fn substitution_score(seq1: *mut c_char, seq2: *mut c_char) -> c_double;
    pub fn substitution_hosei(seq1: *mut c_char, seq2: *mut c_char) -> c_double;
    pub fn substitution_nid(seq1: *mut c_char, seq2: *mut c_char) -> c_double;
    pub fn naivepairscore11(seq1: *mut c_char, seq2: *mut c_char, penal: c_int) -> c_double;
    pub fn naivepairscorefast(
        seq1: *mut c_char,
        seq2: *mut c_char,
        skip1: *mut c_int,
        skip2: *mut c_int,
        penal: c_int,
    ) -> c_double;
    pub fn makeskiptable(n: c_int, skip: *mut *mut c_int, seq: *mut *mut c_char);

    // -- FFT --
    pub fn fft(n: c_int, x: *mut Fukusosuu, dum: c_int) -> c_int;

    // -- Tree construction --
    pub fn nj(nseq: c_int, omtx: *mut *mut c_double, topol: *mut *mut *mut c_int, dis: *mut *mut c_double);
    pub fn upg2(nseq: c_int, eff: *mut *mut c_double, topol: *mut *mut *mut c_int, len: *mut *mut c_double);
    pub fn veryfastsupg_int(nseq: c_int, oeff: *mut *mut c_int, topol: *mut *mut *mut c_int, len: *mut *mut c_double);
    pub fn counteff_simple_double(nseq: c_int, topol: *mut *mut *mut c_int, len: *mut *mut c_double, node: *mut c_double);

    // -- Pairwise alignment --
    pub fn G__align11(
        scoringmtx: *mut *mut c_double,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        alloclen: c_int,
        headgp: c_int,
        tailgp: c_int,
    ) -> c_double;
    pub fn L__align11(
        scoringmtx: *mut *mut c_double,
        scoreoffset: c_double,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        alloclen: c_int,
        off1pt: *mut c_int,
        off2pt: *mut c_int,
    ) -> c_double;

    // -- Tree weighting --
    pub fn treeCnv(
        stopol: *mut Node,
        locnseq: c_int,
        topol: *mut *mut *mut c_int,
        len: *mut *mut c_double,
        bw: *mut *mut c_double,
    );
    pub fn calcBranchWeight(
        bw: *mut *mut c_double,
        locnseq: c_int,
        stopol: *mut Node,
        topol: *mut *mut *mut c_int,
        len: *mut *mut c_double,
    );
    pub fn weightFromABranch(
        nseq: c_int,
        result: *mut c_double,
        stopol: *mut Node,
        topol: *mut *mut *mut c_int,
        step: c_int,
        LorR: c_int,
    );
    pub fn fastconjuction_noname(
        memlist: *mut c_int,
        seq: *mut *mut c_char,
        aseq: *mut *mut c_char,
        peff: *mut c_double,
        eff: *mut c_double,
        d: *mut c_char,
        mineff: c_double,
        oritotal: *mut c_double,
    ) -> c_int;
    pub fn cpmx_calc_new(
        seq: *mut *mut c_char,
        cpmx: *mut *mut c_double,
        eff: *mut c_double,
        lgth: c_int,
        clus: c_int,
    );
    pub fn st_OpeningGapCount(
        ogcp: *mut c_double,
        clus: c_int,
        seq: *mut *mut c_char,
        eff: *mut c_double,
        len: c_int,
    );
    pub fn st_FinalGapCount(
        fgcp: *mut c_double,
        clus: c_int,
        seq: *mut *mut c_char,
        eff: *mut c_double,
        len: c_int,
    );
    pub fn gapcountf(
        freq: *mut c_double,
        seq: *mut *mut c_char,
        nseq: c_int,
        eff: *mut c_double,
        lgth: c_int,
    );
    pub fn MSalignmm(
        n_dynamicmtx: *mut *mut c_double,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        icyc: c_int,
        jcyc: c_int,
        alloclen: c_int,
        sgap1: *mut c_char,
        sgap2: *mut c_char,
        egap1: *mut c_char,
        egap2: *mut c_char,
        chudanpt: *mut c_int,
        chudanref: c_int,
        chudanres: *mut c_int,
        headgp: c_int,
        tailgp: c_int,
        cpmxchild0: *mut *mut *mut c_double,
        cpmxchild1: *mut *mut *mut c_double,
        cpmxresult: *mut *mut *mut c_double,
        orieff1: c_double,
        orieff2: c_double,
    ) -> c_double;
    pub fn A__align(
        n_dynamicmtx: *mut *mut c_double,
        penalty: c_int,
        penalty_ex: c_int,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        icyc: c_int,
        jcyc: c_int,
        alloclen: c_int,
        constraint: c_int,
        impmatch: *mut c_double,
        sgap1: *mut c_char,
        sgap2: *mut c_char,
        egap1: *mut c_char,
        egap2: *mut c_char,
        chudanpt: *mut c_int,
        chudanref: c_int,
        chudanres: *mut c_int,
        headgp: c_int,
        tailgp: c_int,
        firstmem: c_int,
        calledbyfulltreebase: c_int,
        cpmxchild0: *mut *mut *mut c_double,
        cpmxchild1: *mut *mut *mut c_double,
        cpmxresult: *mut *mut *mut c_double,
        orieff1: c_double,
        orieff2: c_double,
    ) -> c_double;
    pub fn imp_match_init_strict(
        imp: *mut c_double,
        clus1: c_int,
        clus2: c_int,
        lgth1: c_int,
        lgth2: c_int,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        eff1_kozo: *mut c_double,
        eff2_kozo: *mut c_double,
        localhom: *mut *mut *mut LocalHom,
        swaplist: *mut c_char,
        forscore: c_int,
        orinum1: *mut c_int,
        orinum2: *mut c_int,
        uselh: *mut c_int,
        seedinlh1: *mut c_int,
        seedinlh2: *mut c_int,
        nodeid: c_int,
        nfiles: c_int,
    );
    pub fn imp_match_out_sc(i1: c_int, j1: c_int) -> c_double;
    pub fn intergroup_score(
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        clus1: c_int,
        clus2: c_int,
        len: c_int,
        value: *mut c_double,
    );

    // -- FFT / segment detection --
    pub fn alignableReagion(
        clus1: c_int,
        clus2: c_int,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        seg: *mut Segment,
    ) -> c_int;

    pub fn searchAnchors(
        nseq: c_int,
        seq: *mut *mut c_char,
        seg: *mut Segment,
    ) -> c_int;

    pub fn fixed_musclesupg_double_treeout(
        nseq: c_int,
        eff: *mut *mut c_double,
        topol: *mut *mut *mut c_int,
        len: *mut *mut c_double,
        name: *mut *mut c_char,
    );

    pub fn Falign(
        whichmtx: *mut *mut c_int,
        scoringmatrices: *mut *mut *mut c_double,
        n_dynamicmtx: *mut *mut c_double,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        eff1s: *mut *mut c_double,
        eff2s: *mut *mut c_double,
        clus1: c_int,
        clus2: c_int,
        alloclen: c_int,
        fftlog: *mut c_int,
        chudanpt: *mut c_int,
        chudanref: c_int,
        chudanres: *mut c_int,
    ) -> c_double;

    // -- Utility --
    pub fn seqlen(seq: *mut c_char) -> c_int;
    pub fn commongappick(nseq: c_int, seq: *mut *mut c_char);
    pub fn reporterr(str: *const c_char, ...);
    pub fn ErrorExit(message: *mut c_char);
}

// ---------------------------------------------------------------------------
// Global variable accessors
// ---------------------------------------------------------------------------
//
// MAFFT's C code uses ~150 extern globals. Rather than declare every one,
// we expose the most critical ones here. Additional globals can be added
// as migration progresses.

unsafe extern "C" {
    pub static mut njob: c_int;
    pub static mut nlenmax: c_int;
    pub static mut dorp: c_int;
    pub static mut alg: c_char;
    pub static mut nthread: c_int;
    pub static mut niter: c_int;
    pub static mut weight: c_int;
    pub static mut utree: c_int;
    pub static mut constraint: c_int;
    pub static mut nblosum: c_int;
    pub static mut TMorJTT: c_int;
    pub static mut fmodel: c_int;
    pub static mut scoremtx: c_int;

    pub static mut penalty: c_int;
    pub static mut ppenalty: c_int;
    pub static mut penalty_dist: c_int;
    pub static mut ppenalty_dist: c_int;
    pub static mut penalty_ex: c_int;
    pub static mut ppenalty_ex: c_int;
    pub static mut penalty_OP: c_int;
    pub static mut ppenalty_OP: c_int;
    pub static mut penalty_EX: c_int;
    pub static mut ppenalty_EX: c_int;
    pub static mut offset: c_int;
    pub static mut poffset: c_int;
    pub static mut offsetFFT: c_int;
    pub static mut offsetLN: c_int;
    pub static mut penaltyLN: c_int;
    pub static mut penalty_exLN: c_int;
    pub static mut RNAppenalty: c_int;
    pub static mut RNAppenalty_ex: c_int;
    pub static mut RNApthr: c_int;
    pub static mut nevermemsave: c_int;
    pub static mut disp: c_int;

    pub static mut use_fft: c_char;
    pub static mut force_fft: c_char;
    pub static mut fftscore: c_int;
    pub static mut fftWinSize: c_int;
    pub static mut fftThreshold: c_int;
    pub static mut fftRepeatStop: c_int;
    pub static mut fftNoAnchStop: c_int;
    pub static mut fftkeika: c_int;
    pub static mut kobetsubunkatsu: c_int;
    pub static mut divWinSize: c_int;
    pub static mut divThreshold: c_int;
    pub static mut score_check: c_int;
    pub static mut bunkatsu: c_int;
    pub static mut fastathreshold: c_double;
    pub static mut penalty_shift_factor: c_double;
    pub static mut cooling: c_int;
    pub static mut legacygapcost: c_int;
    pub static mut consweight_multi: c_double;
    pub static mut specificityconsideration: c_double;
    pub static mut maxdistclass: c_int;
    pub static mut randomseed: c_int;
    pub static mut scmtd: c_int;
    pub static mut refine: c_int;
    pub static mut check: c_int;
    pub static mut cut: c_double;
    pub static mut intop: c_int;
    pub static mut intree: c_int;
    pub static mut devide: c_int;
    pub static mut rnakozo: c_int;
    pub static mut rnaprediction: c_char;

    pub static mut outgap: c_int;
    pub static mut kimuraR: c_int;
    pub static mut pamN: c_int;
    pub static mut treemethod: c_int;
    pub static mut sueff_global: c_double;
    pub static mut minimumweight: c_double;

    pub static mut amino_dis: *mut *mut c_int;
    pub static mut n_dis: *mut *mut c_int;
    pub static mut n_disFFT: *mut *mut c_int;
    pub static mut amino_dis_consweight_multi: *mut *mut c_double;
    pub static mut n_dis_consweight_multi: *mut *mut c_double;
    pub static mut amino: [u8; 0x100];
    pub static mut polarity: [c_double; 0x100];
    pub static mut volume: [c_double; 0x100];
    pub static mut ribosumdis: [[c_int; 37]; 37];
    pub static mut amino_n: [c_int; 0x100];
    pub static mut amino_grp: [c_char; 0x100];
    pub static mut nalphabets: c_int;
    pub static mut nscoredalphabets: c_int;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr::addr_of;

    #[test]
    fn test_linking_and_init() {
        // Verify that we can call into the C library.
        // initglobalvariables() sets defaults for all C globals.
        unsafe {
            initglobalvariables();
            assert_eq!(addr_of!(nthread).read(), 1);
            assert_eq!(addr_of!(outgap).read(), 1);
        }
    }
}
