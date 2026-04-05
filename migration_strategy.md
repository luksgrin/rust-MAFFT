# MAFFT C-to-Rust Migration Strategy

## Overview

This document outlines an incremental migration plan for rewriting the MAFFT multiple
sequence alignment tool from C to Rust. The strategy prioritizes correctness at every step:
each phase produces a hybrid C/Rust binary that passes the existing test suite before
proceeding.

The core idea is **bottom-up replacement via C-FFI boundaries**. We start with leaf modules
(no internal dependents), expose them to C via `extern "C"` functions, swap the linker
inputs, and work our way up until the orchestration layer is the last piece rewritten.

---

## Current Architecture (Summary)

```
mafft (shell script orchestrator)
  |
  +-- disttbfast   (distance matrix + guide tree + progressive alignment)
  +-- tbfast       (tree-based fast alignment + optional iterative refinement)
  +-- dvtditr      (tree-dependent iterative refinement)
  +-- splittbfast  (memory-efficient splitting variant)
  +-- addsingle    (add sequences to existing alignment)
  +-- dndfast7     (distance-only calculation)
  +-- f2cl, seq2regtable, regtable2seq, ...  (utilities)
```

All binaries link against a common set of ~20 C object files. The central header
`mltaln.h` defines shared data structures and ~400 extern globals.

**Key coupling points:**
- Global scoring matrices initialized by `constants()` at startup
- `LocalHom` linked lists threaded through every alignment function
- Manual 2D/3D matrix allocation (`AllocateXxxMtx` / `FreeXxxMtx`)
- Thread-local storage via `TLS` macro + pthreads
- Recursive tree traversal with raw pointer arithmetic

---

## Guiding Principles

1. **Always shippable.** After every phase the project compiles, links, and passes tests.
   Never break the existing C pipeline while Rust code is incomplete.

2. **FFI as the seam.** Each migrated module exposes `extern "C"` functions matching the
   original C signatures. The C side includes a thin header that calls into Rust. This
   lets us swap one `.o` at a time in the Makefile.

3. **Test-driven confidence.** Before rewriting a module, write Rust integration tests
   that call the *C* version through FFI and capture reference outputs. The Rust
   replacement must produce bit-identical results.

4. **Eliminate globals gradually.** Each Rust module owns its state in structs. At FFI
   boundaries we translate between the C global world and Rust's owned-data world. By the
   final phase, globals are gone.

5. **No premature optimization.** First make it correct. Rust's `Vec`, `ndarray`, and
   `rayon` will naturally replace manual allocation and pthreads, but only after
   correctness is proven.

---

## Project Setup (Phase 0)

**Goal:** Establish the hybrid build infrastructure.

### Steps

1. Initialize a Cargo workspace at the repo root:
   ```
   rust-MAFFT/
     Cargo.toml          (workspace)
     core/               (existing C code, untouched)
     crates/
       mafft-sys/        (raw FFI bindings to remaining C code)
       mafft-types/      (shared Rust types: LocalHom, Node, Segment, etc.)
       mafft-io/         (Phase 1)
       mafft-scoring/    (Phase 2)
       mafft-fft/        (Phase 3)
       mafft-align/      (Phase 4-5)
       mafft-tree/       (Phase 4)
       mafft-core/       (Phase 6)
       mafft-bin/        (Phase 7 -- final binaries)
     scripts/            (existing shell wrappers, adapted last)
     test/               (existing test data)
   ```

2. Create `mafft-sys` using a `build.rs` that compiles the C sources with the `cc` crate.
   This gives us a working baseline: Rust tests can call any C function.

3. Create `mafft-types` with pure-Rust equivalents of the core structs from `mltaln.h`:
   - `LocalHom` -> `LocalHomology { regions: Vec<HomologyRegion> }` (no linked list)
   - `Node` -> enum-based tree with arena allocation
   - `Segment` -> `AlignmentSegment` with index-based pairing
   - `Fukusosuu` -> `num::Complex<f64>` (or a minimal custom type)
   - Scoring matrices as owned 2D arrays
   - Conversion traits (`From<*const C_LocalHom>`, etc.) for FFI boundaries

4. Set up CI: `cargo test` + existing C test suite (run the shell-script tests that
   compare output against `test/` references).

**Deliverable:** `cargo test` passes. All existing C binaries still build and produce
correct output. No C code modified yet.

---

## Phase 1 -- I/O (`mafft-io`)

**Target files:** `io.c` (6.5K lines), parts of `mtxutl.c`

**Why first:** I/O is a leaf dependency -- many modules call it, but it calls almost
nothing else. It has a clean interface (`PreRead`, `Read`, `readData`, `Write`). Replacing
it also lets us define Rust-native sequence types early.

### Rust API

```rust
pub struct SequenceSet {
    pub names: Vec<String>,
    pub sequences: Vec<Vec<u8>>,  // raw residue codes
    pub comments: Vec<String>,
}

pub fn read_fasta(path: &Path) -> Result<SequenceSet>;
pub fn write_fasta(seqs: &SequenceSet, path: &Path, format: OutputFormat) -> Result<()>;
```

### Steps

1. Write `mafft-io` with FASTA parser and writer.
2. Write integration tests that round-trip every file in `test/`.
3. Expose `extern "C"` wrappers: `PreRead_rs`, `Read_rs`, `writeData_rs`.
4. In `core/Makefile`, replace `io.o` with the Rust static library.
5. Run full test suite.

### Validation

- Byte-for-byte output comparison on `test/` inputs.
- Edge cases: empty sequences, very long names, non-standard characters.

---

## Phase 2 -- Scoring Matrices & Constants (`mafft-scoring`)

**Target files:** `constants.c`, `defs.c`, parts of `mltaln9.c` (matrix init)

**Why now:** Every algorithm depends on scoring matrices. Migrating them early means
subsequent Rust modules can use native Rust types instead of reading C globals.

### Rust API

```rust
pub struct ScoringContext {
    pub amino_dis: Array2<i32>,        // amino acid substitution matrix
    pub n_dis: Array2<i32>,            // nucleotide substitution matrix
    pub n_dis_fft: Array2<i32>,        // nucleotide matrix for FFT
    pub polarity: [f64; 256],
    pub volume: [f64; 256],
    pub amino_map: [u8; 256],          // char -> internal index
    pub penalty: GapPenalties,
    pub seq_type: SeqType,             // DNA or Protein
}

pub fn load_scoring(model: Model, seq_type: SeqType) -> ScoringContext;
```

### Steps

1. Port BLOSUM, JTT, TM, and DNA matrices as const arrays in Rust.
2. Port the `constants()` initialization logic.
3. Expose via FFI: a single `init_scoring_context()` that populates C-visible globals
   *from* the Rust-owned data (backward compat while C code still reads globals).
4. Replace `constants.o` and `defs.o` in the Makefile.
5. Run full test suite.

---

## Phase 3 -- FFT Engine (`mafft-fft`)

**Target files:** `fft.c`, `fftFunctions.c`

**Why now:** The FFT module is mathematically self-contained (complex arithmetic +
Cooley-Tukey). It depends only on scoring matrices (Phase 2). It is the signature
algorithmic innovation of MAFFT and benefits greatly from Rust's safety guarantees around
buffer indexing.

### Rust API

```rust
pub fn fft_homology_search(
    seq1: &[u8],
    seq2: &[u8],
    scoring: &ScoringContext,
    params: &FftParams,
) -> Vec<HomologyRegion>;
```

### Steps

1. Port the Cooley-Tukey FFT (`fft.c`).
2. Port the sequence-to-vector conversion and correlation scoring (`fftFunctions.c`).
3. Validate against C output on all test pairs (exact floating-point match).
4. Swap `fft.o` and `fftFunctions.o` in the Makefile.

### Notes

- Consider using `rustfft` crate for the core transform -- it's well-tested and SIMD-
  optimized. Validate that it produces identical results to the hand-rolled C FFT.
- The `Fukusosuu` (complex number) type maps directly to `num::Complex<f64>`.

---

## Phase 4 -- Pairwise Alignment & Tree Building (`mafft-align`, `mafft-tree`)

**Target files:**
- Pairwise: `Galign11.c`, `Lalign11.c`, `genalign11.c`, `Salignmm.c`, `Dalignmm.c`,
  `Lalignmm.c`, `partSalignmm.c`
- Tree: `tddis.c`, NJ/UPGMA portions of `mltaln9.c`

These can be developed **in parallel** by two contributors since they have minimal
cross-dependency.

### Phase 4a -- Pairwise Alignment (`mafft-align`)

```rust
pub fn global_align(s1: &[u8], s2: &[u8], scoring: &ScoringContext) -> Alignment;
pub fn local_align(s1: &[u8], s2: &[u8], scoring: &ScoringContext) -> Vec<LocalHit>;
pub fn semiglobal_align(...) -> Alignment;

// Group-to-group (profile) alignment
pub fn profile_align(
    group1: &AlignedGroup,
    group2: &AlignedGroup,
    scoring: &ScoringContext,
) -> Alignment;
```

#### Steps

1. Port Needleman-Wunsch (`Galign11.c`) with affine gaps.
2. Port Smith-Waterman (`Lalign11.c`, `Lalignmm.c`).
3. Port DNA-specific alignment (`Dalignmm.c`).
4. Port generalized affine gap model (`genalign11.c`) for E-INS-i.
5. Port group-to-group alignment (`Salignmm.c`, `partSalignmm.c`).
6. For each: validate against C on test inputs, then swap `.o` files.

#### Notes

- Replace the manual 2D DP matrix allocation (`AllocateFloatMtx`) with `Vec<Vec<f32>>`
  or a flat `Vec<f32>` with stride-based indexing.
- Replace thread-local `commonIP`/`commonJP` traceback arrays with function-local
  allocations (Rust's allocator is fast enough; profile before optimizing).

### Phase 4b -- Distance & Tree Construction (`mafft-tree`)

```rust
pub fn compute_distance_matrix(
    seqs: &SequenceSet,
    method: DistanceMethod,  // KTuple, Local, Global
    scoring: &ScoringContext,
) -> DistanceMatrix;

pub fn neighbor_joining(dist: &DistanceMatrix) -> GuideTree;
pub fn upgma(dist: &DistanceMatrix) -> GuideTree;
pub fn parttree(seqs: &SequenceSet, ...) -> GuideTree;  // for 10K+ seqs
```

#### Steps

1. Port k-tuple distance calculation.
2. Port neighbor-joining (`nj()` in `mltaln9.c`).
3. Port UPGMA variants (`upg2()`, `veryfastsupg_int()`).
4. Port PartTree for large datasets.
5. Validate tree topologies against C output.

---

## Phase 5 -- Multi-Sequence & FFT-Accelerated Alignment

**Target files:** `MSalignmm.c`, `MSalign11.c`, `Falign.c`, `Falign_localhom.c`,
`SAalignmm.c`, `pairlocalalign.c`, `rna.c`

**Why after Phase 4:** These depend on pairwise alignment (Phase 4a), FFT (Phase 3), and
tree structures (Phase 4b).

### Steps

1. Port `MSalignmm` -- progressive multi-sequence alignment (the workhorse).
2. Port `Falign` -- FFT-accelerated alignment using Phase 3's FFT engine.
3. Port `Falign_localhom` -- FFT alignment with local homology constraints.
4. Port `pairlocalalign` -- the constraint-building pairwise step for L-INS-i / E-INS-i.
5. Port `SAalignmm` -- structure-aware alignment.
6. Port `rna.c` -- RNA-specific scoring adjustments.

### Validation

At this point we can run end-to-end alignment on test sequences through the Rust
implementations of every algorithmic layer. Compare output alignments character-by-character
against the C reference.

---

## Phase 6 -- Iterative Refinement & Orchestration Core (`mafft-core`)

**Target files:** `tditeration.c`, `addfunctions.c`, large portions of `mltaln9.c`

This is the most complex phase. `TreeDependentIteration` is the heart of MAFFT's accuracy
modes (L-INS-i, G-INS-i, E-INS-i). It orchestrates repeated re-alignment guided by the
tree and local homology constraints.

### Rust API

```rust
pub struct MafftEngine {
    scoring: ScoringContext,
    params: AlignmentParams,
}

impl MafftEngine {
    pub fn progressive_align(
        &self,
        seqs: &SequenceSet,
        tree: &GuideTree,
    ) -> MultipleAlignment;

    pub fn iterative_refine(
        &self,
        alignment: &mut MultipleAlignment,
        tree: &GuideTree,
        local_hom: &LocalHomologyTable,
        max_iterations: usize,
    ) -> RefinementResult;
}
```

### Steps

1. Port `TreeDependentIteration` loop logic.
2. Port weight calculation and score convergence checks.
3. Port `addfunctions.c` (--add / --addfragments functionality).
4. This phase eliminates the last dependency on `mltaln9.o`.

---

## Phase 7 -- Entry Points & Binary Consolidation (`mafft-bin`)

**Target files:** `dvtditr.c`, `disttbfast.c`, `tbfast.c`, `splittbfast.c`, `addsingle.c`

### Steps

1. Port command-line argument parsing. Use `clap` for a proper CLI.
2. Port each entry point as a Rust binary that calls into `mafft-core`.
3. Consider **consolidating** binaries: instead of 5+ separate executables, produce a
   single `mafft` binary with subcommands, simplifying distribution.
4. Rewrite the shell orchestrator (`mafft.tmpl`) in Rust as the top-level CLI, eliminating
   the shell dependency.
5. Port small utilities (`f2cl`, `seq2regtable`, `regtable2seq`, etc.).

### Deliverable

A single `mafft` binary (or a small set) written entirely in Rust, with no C dependencies.

---

## Phase 8 -- Cleanup & Optimization

1. **Remove `mafft-sys`** -- no more C code to bind to.
2. **Parallelism:** Replace pthreads with `rayon` for data-parallel distance matrix and
   pairwise alignment computation.
3. **SIMD:** Use `std::simd` (or `pulp`/`simdeez`) for inner DP loops and FFT.
4. **Memory:** Profile and optimize. Rust's ownership model should reduce peak memory vs.
   the C version's conservative allocation patterns.
5. **API crate:** Expose `mafft-core` as a library crate on crates.io so other Rust
   bioinformatics tools can call MAFFT programmatically.
6. **Python bindings:** Use PyO3 to create `pymafft` for the Python bioinformatics
   ecosystem.

---

## Phase Summary

| Phase | Module(s) | C Files Replaced | Risk | Est. Complexity |
|-------|-----------|-------------------|------|-----------------|
| 0 | Build infra, types | none | Low | Medium |
| 1 | I/O | io.c, mtxutl.c | Low | Low |
| 2 | Scoring matrices | constants.c, defs.c | Low | Low |
| 3 | FFT engine | fft.c, fftFunctions.c | Medium | Medium |
| 4a | Pairwise alignment | Galign11, Lalign11, Salignmm, Dalignmm, etc. | Medium | High |
| 4b | Tree building | tddis.c, NJ/UPGMA in mltaln9 | Medium | Medium |
| 5 | Multi-seq alignment | MSalignmm, Falign, pairlocalalign, rna | High | High |
| 6 | Iterative refinement | tditeration, addfunctions, mltaln9 | High | High |
| 7 | Entry points & CLI | dvtditr, disttbfast, tbfast, etc. | Medium | Medium |
| 8 | Cleanup & optimization | (remove all C) | Low | Medium |

---

## Risk Mitigation

**Numerical divergence.** Floating-point results may differ between C and Rust due to
compiler optimizations and operation ordering. Mitigation: use `#[cfg(test)]` comparison
with epsilon tolerance; for critical paths, match C's exact operation order initially.

**Global state.** The C code uses ~400 extern globals. Mitigation: Phase 0's FFI layer
reads/writes globals through accessor functions. Each Rust module uses owned structs
internally. Globals shrink as modules are ported.

**`mltaln9.c` is 16K lines.** This monolith contains tree iteration, NJ, scoring, memory
helpers, and more. Mitigation: do *not* port it as one unit. Extract functions into the
appropriate Rust crate as they're needed (tree code to `mafft-tree`, scoring to
`mafft-align`, etc.). What remains gets smaller each phase.

**Test coverage.** The existing test suite is small (a few reference alignments).
Mitigation: before each phase, generate additional reference outputs by running the C
version on diverse inputs. Store these as golden files for Rust tests.

---

## Dependency Graph

```
mafft-bin (Phase 7)
  |
  +-- mafft-core (Phase 6)
        |
        +-- mafft-align (Phase 4a, 5)
        |     |
        |     +-- mafft-fft (Phase 3)
        |     |     |
        |     |     +-- mafft-scoring (Phase 2)
        |     |
        |     +-- mafft-scoring (Phase 2)
        |
        +-- mafft-tree (Phase 4b)
        |     |
        |     +-- mafft-scoring (Phase 2)
        |
        +-- mafft-io (Phase 1)
        |
        +-- mafft-types (Phase 0)

mafft-sys (Phase 0, removed in Phase 8)
  |
  +-- core/*.c (shrinks each phase)
```
