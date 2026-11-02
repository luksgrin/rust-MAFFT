# `--add` follow-up references

Stdout (and, where noted, the `<addfile>.map` file) of the in-tree C
build (`make -C mafft-upstream/core`, arm64 / Apple clang), run with
`LC_ALL=C --thread 1`. Each input is an aligned existing set plus
unaligned sequences to add:

* `uniprot.*` — 33 aligned rows + 3 new, all with UniProt-style
  `sp|…|…` / `tr|…|…` titles and a non-standard residue (`#`) in several
  rows, so only `--anysymbol` accepts them.
* `sample33.*` — `sample.first36` split 33 + 3; the existing rows keep
  all-gap columns.
* `bb12019.*` — BAliBASE BB12019 split 3 + 2.
* `bb12019g.*`, `sample33g.*`, `sample20g.*`, `bb20027g.*` — BB12019
  3 + 2, `sample.first36` 33 + 3 and 20 + 16, BAliBASE BB20027 25 + 4,
  with the existing alignment's all-gap columns removed. Default
  `--add` runs two guide-tree cycles (`disttbfast -E 2`); these pin
  that, the 16-sequence add (merges among added sequences) and the
  linkage options.

```bash
make -C mafft-upstream/core
export MAFFT_BINARIES=$PWD/mafft-upstream/binaries LC_ALL=C
M=mafft-upstream/scripts/mafft
O=crates/mafft-bin/tests/fixtures/add_followups
T=$(mktemp -d)
# run <tag> <case> <flag> [extra...]: copy the inputs aside, because
# --mapout / --compactmapout write <addfile>.map next to the added file.
run() {
  tag=$1 c=$2 flag=$3; shift 3
  cp $O/$c.existing.fa $O/$c.new.fa $T/
  $M --quiet --thread 1 $flag $T/$c.new.fa "$@" $T/$c.existing.fa > $O/$tag.aln
  if [ -f $T/$c.new.fa.map ]; then mv $T/$c.new.fa.map $O/$tag.map; fi
}
run uniprot.add.anysymbol.keeplength            uniprot   --add          --anysymbol --keeplength
run uniprot.add.anysymbol.compactmapout         uniprot   --add          --anysymbol --compactmapout
run uniprot.addfragments.anysymbol.keeplength   uniprot   --addfragments --anysymbol --keeplength
run sample33.add.keeplength                     sample33  --add          --keeplength
run sample33.add.mapout                         sample33  --add          --mapout
run sample33.addfragments.compactmapout         sample33  --addfragments --compactmapout
run sample33.addfragments.reorder               sample33  --addfragments --reorder
run bb12019.addfragments.reorder                bb12019   --addfragments --reorder
run bb12019.addfragments.averagelinkage         bb12019   --addfragments --averagelinkage
run bb12019.addfragments.minimumlinkage         bb12019   --addfragments --minimumlinkage
run bb12019g.add                                bb12019g  --add
run bb12019g.add.keeplength                     bb12019g  --add          --keeplength
run bb12019g.add.mapout                         bb12019g  --add          --mapout
run bb12019g.add.reorder                        bb12019g  --add          --reorder
run bb12019g.add.retree1                        bb12019g  --add          --retree 1
run sample33g.add                               sample33g --add
run sample33g.add.compactmapout                 sample33g --add          --compactmapout
run sample33g.add.reorder                       sample33g --add          --reorder
run sample20g.add                               sample20g --add
run sample20g.add.retree1                       sample20g --add          --retree 1
run sample20g.add.keeplength                    sample20g --add          --keeplength
run bb20027g.add                                bb20027g  --add
run bb20027g.add.averagelinkage                 bb20027g  --add          --averagelinkage
run bb20027g.add.mapout                         bb20027g  --add          --mapout
```
