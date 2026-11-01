# `--add` variant references (issue #50)

Stdout of the in-tree C build (`make -C mafft-upstream/core`, default
flags, arm64 / Apple clang) for five `--add` variants. Two inputs:

* `sample.add6.*` — the canonical 30 + 6 fixture
  (`crates/mafft-core/tests/fixtures/sample.first30.fftns2.aln` plus
  `sample.last6_for_add.fa`).
* `sentinel.*` — the issue #40 reproducer
  (`crates/mafft-bin/tests/fixtures/memsafety/add_sentinel_existing.fa`
  plus `add_sentinel_new.fa`, 30 aligned rows + 6 new).

`--addfragments` is the script's default `distance=multi` path
(`dndpre -M 2`, `pairlocalalign -Y`, `addsingle -d -l 2.7`), not the
large-input 6-mer fallback and not `--addfragments --localpair` (C
reads off the end of a buffer on some inputs of that combination; see
`docs/architecture/memory-safety.md`).
The added file is the same six full-length sequences used for `--add`.

```bash
make -C mafft-upstream/core
export MAFFT_BINARIES=$PWD/mafft-upstream/binaries
M=mafft-upstream/scripts/mafft
F=crates/mafft-core/tests/fixtures
S=crates/mafft-bin/tests/fixtures/memsafety
O=crates/mafft-bin/tests/fixtures/add_variants
# Copy the inputs aside first: --compactmapout writes <addfile>.map
# next to the added sequences.
$M --quiet --add $F/sample.last6_for_add.fa --maxiterate 2 $F/sample.first30.fftns2.aln > $O/sample.add6.maxiterate2.aln
$M --quiet --add $F/sample.last6_for_add.fa --reorder $F/sample.first30.fftns2.aln > $O/sample.add6.reorder.aln
$M --quiet --addfragments $F/sample.last6_for_add.fa $F/sample.first30.fftns2.aln > $O/sample.add6.addfragments.aln
$M --quiet --addfragments $F/sample.last6_for_add.fa --keeplength $F/sample.first30.fftns2.aln > $O/sample.add6.addfragments.keeplength.aln
$M --quiet --add $F/sample.last6_for_add.fa --compactmapout $F/sample.first30.fftns2.aln > $O/sample.add6.compactmapout.aln
$M --quiet --add $S/add_sentinel_new.fa --maxiterate 2 $S/add_sentinel_existing.fa > $O/sentinel.maxiterate2.aln
$M --quiet --add $S/add_sentinel_new.fa --reorder $S/add_sentinel_existing.fa > $O/sentinel.reorder.aln
$M --quiet --addfragments $S/add_sentinel_new.fa $S/add_sentinel_existing.fa > $O/sentinel.addfragments.aln
$M --quiet --addfragments $S/add_sentinel_new.fa --keeplength $S/add_sentinel_existing.fa > $O/sentinel.addfragments.keeplength.aln
$M --quiet --add $S/add_sentinel_new.fa --compactmapout $S/add_sentinel_existing.fa > $O/sentinel.compactmapout.aln
```
