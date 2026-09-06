# Panaroo tiny/core parity corpus — DNA gene clusters

Five real per-gene nucleotide clusters (4 sequences each, ~300–2300 bp) from
the panaroo-rs tiny/core parity corpus, contributed via issue #1 by Johan
Henriksson (mahogny/rust-MAFFT, commit 0b3955d). They are the exact workload
panaroo-rs hands to MAFFT once per gene family:

```bash
mafft --auto --adjustdirection --thread 1 --nuc <cluster.fa>
```

The `.fa` inputs are taken verbatim from that contribution. The expected
outputs were **not** taken from it: every `.expected` here was regenerated on
2026-09-06 from the in-tree reference build (the pinned `mafft-upstream`
source compiled with `make -C mafft-upstream/core`, upstream Makefile default
flags, arm64/clang — see `docs/architecture/byte-identity.md`) and then found
byte-identical to the fork's copies, which were captured on x86-64 without FMA.
These inputs are therefore not sensitive to the floating-point contraction
policy, so a single variant serves both platforms; if that ever changes, add
the fixture to `crates/mafft-core/tests/fixtures/policy_fixtures.tsv` and
regenerate with `scripts/regen_policy_fixtures.sh`.

Two reference files per cluster:

- `<name>.expected` — `mafft --quiet --auto --adjustdirection --thread 1 --nuc`
  (the panaroo-rs invocation);
- `<name>.nothread.expected` — same without `--thread 1`, i.e. C's
  single-threaded refinement path. Currently byte-identical to `.expected`
  for all five clusters; kept separately because the two C code paths differ
  in convergence semantics and could diverge on other inputs.

Regenerate with:

```bash
make -C mafft-upstream/core
export MAFFT_BINARIES=$PWD/mafft-upstream/binaries
D=crates/mafft-bin/tests/fixtures/panaroo_tiny_dna_clusters
for f in "$D"/*.fa; do
  mafft-upstream/scripts/mafft --quiet --auto --adjustdirection --thread 1 --nuc "$f" > "${f%.fa}.expected"
  mafft-upstream/scripts/mafft --quiet --auto --adjustdirection --nuc "$f" > "${f%.fa}.nothread.expected"
done
```

Consumed by the `panaroo_*` tests in `crates/mafft-bin/tests/c_parity_dna.rs`.
