# Known Issues and Future Work

## Per-group gap stripping in progressive alignment

**Status**: Not implemented (using global-all-gap stripping instead)  
**Priority**: Medium (performance optimization)  
**Location**: `crates/mafft-core/src/progressive.rs`, `merge_step()`

### Background

The C code (`commongappick()` in `disttbfast.c`) strips columns that are all-gap within each group **independently** before profile alignment. This means group1's profile only contains columns where at least one group1 member has a residue, and likewise for group2. The two stripped profiles may have different lengths.

Our Rust code currently strips only columns that are all-gap across **all sequences globally**. This is more conservative: it only removes columns where every sequence in the entire MSA has a gap.

### Performance impact

For a 500-column MSA where group1 has 100 residue-containing columns and group2 has 120:
- C's per-group: DP aligns 100 × 120 = 12,000 cells
- Our global: DP aligns ~400 × 400 = 160,000 cells (if 100 columns are globally all-gap)
- Worst case: no globally-all-gap columns → DP aligns 500 × 500 = 250,000 cells

This gap widens with larger alignments and more sequences.

### What was tried

Per-group stripping was implemented and passes unit tests (2-6 sequences) but fails on real data (36 protein sequences from `test/sample`). The two-pass architecture works:
1. Pass 1: Align stripped profiles using `kept1`/`kept2` column index maps.
2. Pass 2: Re-insert stripped columns at original positions.

The bug is in pass 2's interleaving logic when `kept1` and `kept2` have different lengths (which is the whole point of per-group stripping). The re-insertion must correctly interleave:
- Aligned result columns (from ops on stripped profiles)
- Group1-only gap columns (all-gap in group1, content in group2/others)
- Group2-only gap columns (all-gap in group2, content in group1/others)
- Both-gap columns (all-gap in both groups, possible content in others)

The C code handles this in `insertnewgaps()` (`addfunctions.c`, lines 445-673) by walking through aligned results position-by-position, classifying each sequence into group1/group2/other, and copying content or gaps accordingly. The key is that it processes ALL sequences simultaneously, column by column, using the alignment result as the master ordering.

### How to fix

1. Study `insertnewgaps()` in the C code more carefully — specifically how it handles the case where `gap1 != gap2` (different columns stripped from each group).
2. The interleaving needs to merge three column streams:
   - Result columns from the alignment (ordered by ops)
   - Group1-stripped columns (at their original positions relative to kept1)
   - Group2-stripped columns (at their original positions relative to kept2)
3. A correct merge must assign each original column index to exactly one output position, preserving the relative order of all original columns.
4. Test with the 36-sequence `test/sample` dataset, specifically checking sequence 3 (`M92038 chicken green sensitive cone opsin`) which was the first to show corruption.

### Iterative refinement gap bookkeeping

**Status**: Working correctly with the "all sequences split" approach  
**Location**: `crates/mafft-core/src/refinement.rs`

The refinement code splits ALL sequences into two groups at each tree branch (matching C's `OneClusterAndTheOther_fast()`), ensuring no "other" sequences exist during re-alignment. This eliminates the column-space ambiguity problem. The same per-group stripping optimization could be applied here too, but the same interleaving bug would apply.
