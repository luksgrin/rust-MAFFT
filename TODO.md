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

### Root cause (diagnosed)

The cursor model breaks when `kept1 != kept2`. A Match op aligning `kept1[c1]` (original column X) with `kept2[c2]` (original column Y, where X ≠ Y) emits one output column. "Other" sequences follow cursor1 and get content from column X. Their content at column Y is lost because:
- Column Y was consumed by cursor2 (as part of the Match)
- But "other" sequences don't follow cursor2
- Skipping column Y drops "other"'s content there
- Emitting column Y directly would double-emit (it's already part of the Match op for group2)

Multiple approaches were tried:
1. **Two-pass with interleaving** — ops placed at primary column positions, Insert ops interleaved. Failed because Insert columns referencing earlier original positions get placed out of order.
2. **Column classification with consumed flags** — marks all op-consumed columns, skips them in direct emission. Failed because "other" content at consumed columns is lost.
3. **Simple cursor with trailing gap columns** — appends both-gap columns at end. Failed because K1 < width drops "other" content at group1-gap columns.
4. **Conditional Insert: emit "other" content at kept2-only columns** — At Insert positions, "other" gets content from kept2[c2] if that column is NOT in kept1. Failed because the DP can cross-align columns in both kept sets, causing double-emission (357 vs 353 residues).
5. **Template approach (O/=/~ markers)** — Build a template marking each output position as old-column/new-gap. Walk template with shared orig_cursor. Failed because group2 needs a SEPARATE cursor: at Delete positions (O/=), group2 gets gap but orig_cursor advances, causing group2 to skip columns and lose content.
6. **Three-cursor template** — separate g1_cursor and g2_cursor. Difficult because at O/O positions (old column for both), BOTH cursors must advance in sync, but at Delete positions only g1 advances. This causes g2_cursor to fall behind, and the total number of g2 advances (O/O + =/~ = width - Delete + Insert = width - K1 + K2) only equals width when K1 == K2.

### How to fix

The correct approach requires porting C's `insertnewgaps()` from `addfunctions.c` (lines 445-673). The C function uses a fundamentally different model:

1. **It does NOT use alignment ops.** Instead, it compares the aligned result (sequences with '=' markers for new gaps) character by character with the original sequences.
2. **It uses `gaplen[]` and `gapmap[]` arrays** built by `findnewgaps()` and `adjustgapmap()` to determine where new gaps appear in the merged result.
3. **It walks the merged result position by position**, using `posin12` (position in group1/group2's aligned sequence) and `j` (position in group0/other's original sequence) as independent counters.
4. **At new-gap positions**: group0 gets gap characters, group1/group2 advance through their aligned content.
5. **At old-column positions**: all three groups advance through their content.

The key difference from our ops-based model: C's function operates on the RESULT sequences directly (character arrays with embedded gap markers), not on abstract alignment operations. This avoids the cursor-synchronization problem entirely because each group's position counter only advances when that group actually has content at the current position.

To port this:
1. After the profile DP, reconstruct group1 and group2's full-width aligned sequences with '=' markers at new-gap positions.
2. Port `findnewgaps()` to scan the '=' markers and build `gaplen[]`.
3. Port `adjustgapmap()` to build `gapmap[]`.
4. Port the main loop of `insertnewgaps()` using `j` and `posin12` counters.
5. This requires modifying the profile alignment to return actual aligned character sequences (not just ops), or reconstructing them from ops first.

### Iterative refinement gap bookkeeping

**Status**: Working correctly with the "all sequences split" approach  
**Location**: `crates/mafft-core/src/refinement.rs`

The refinement code splits ALL sequences into two groups at each tree branch (matching C's `OneClusterAndTheOther_fast()`), ensuring no "other" sequences exist during re-alignment. This eliminates the column-space ambiguity problem. The same per-group stripping optimization could be applied here too, but the same interleaving bug would apply.
