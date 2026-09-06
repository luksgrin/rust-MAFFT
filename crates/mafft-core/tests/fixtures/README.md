# Test fixtures

Reference data committed in this repo to guard the byte-level parity with
C MAFFT 7.526. Upstream does not ship references for every alignment mode
(notably NW-NS-2 / `--nofft`), so we generate and commit our own.

Regenerating any of these requires the C binaries to be built
(`make -C mafft-upstream/core`).

## Per-policy fixtures: `<name>.fma` / `<name>.nofma`

C MAFFT 7.526 is not bit-reproducible across CPU architectures: the arm64
clang build contracts `a*b+c` into fused multiply-adds (~1030 `fmadd` per
engine binary, see `scripts/fma_census.sh`), the baseline x86-64 gcc build
(bioconda) emits none. Where that changes the output, the fixture is
committed twice and the plain name is removed:

- `<name>.fma` — output of the arm64 C binary (fused).
- `<name>.nofma` — output of the x86-64 C binary (not fused).

`fixture_path_fp(name)` in `end_to_end.rs` returns the variant matching
`mafft_types::fp::CONTRACTS_FMA` (the build's contraction policy, selected by
the `fp-contract-fma` / `fp-contract-none` cargo features, defaulting to
fused on aarch64 only) and falls back to the plain `<name>` when no variant
exists. Do not commit a plain `<name>` alongside a variant pair.

### `sample.bl50.fftns2.fma` / `sample.bl50.fftns2.nofma`

C's alignment output for `mafft --bl 50 --retree 2 --maxiterate 0
mafft-upstream/test/sample`: width 712 from the arm64 macOS binary, width
738 from the x86-64 bioconda binary (the latter contributed by @mahogny,
mahogny/rust-MAFFT). `fftns2_bl50_byte_identical_to_c` asserts against the
variant matching the build's policy.

```bash
mafft --quiet --bl 50 --retree 2 --maxiterate 0 mafft-upstream/test/sample \
  > crates/mafft-core/tests/fixtures/sample.bl50.fftns2.fma    # on arm64
  > crates/mafft-core/tests/fixtures/sample.bl50.fftns2.nofma  # on x86-64
```

## `sample.nwns2`

C's alignment output for `mafft --nofft mafft-upstream/test/sample`. The
`nofft_byte_identical_to_c` test asserts Rust's NW-NS-2 output matches this
byte-for-byte.

```bash
MAFFT_BINARIES=$PWD/mafft-upstream/binaries \
  $PWD/mafft-upstream/scripts/mafft --nofft --quiet $PWD/mafft-upstream/test/sample \
  > crates/mafft-core/tests/fixtures/sample.nwns2
```

## `sample.nwns2.op25`

C's alignment output for `mafft --nofft --op 2.5 mafft-upstream/test/sample`.
The `nofft_op_override_byte_identical_to_c` test asserts Rust's NW-NS-2
output with `--op 2.5` matches this byte-for-byte, guarding the command-line
gap-opening override (`GapParams` application in `crates/mafft-core/src/engine.rs`).

```bash
MAFFT_BINARIES=$PWD/mafft-upstream/binaries \
  $PWD/mafft-upstream/scripts/mafft --nofft --op 2.5 --quiet $PWD/mafft-upstream/test/sample \
  > crates/mafft-core/tests/fixtures/sample.nwns2.op25
```

## `samplerna.nwns2`

C's alignment output for `mafft --nofft mafft-upstream/test/samplerna` (a
5-sequence RNA input). The `rna_nofft_case_insensitive_identical_to_c` test
asserts Rust's output matches this byte-for-byte **after case normalization**
on both sides. C preserves the input lowercase; Rust currently uppercases
before alignment. The alignment columns themselves are identical.

```bash
MAFFT_BINARIES=$PWD/mafft-upstream/binaries \
  $PWD/mafft-upstream/scripts/mafft --nofft --quiet $PWD/mafft-upstream/test/samplerna \
  > crates/mafft-core/tests/fixtures/samplerna.nwns2
```

## `sample.nwns2.steps`

C's per-step progressive-alignment scores for `--nofft` on the same input
(one `RDBG step clus1 clus2 width score` line per merge — 70 lines covering
both retree passes). The `nofft_per_step_matches_c` test asserts Rust's
per-step RDBG trace matches this line-by-line.

This is a finer-grained regression guard than `sample.nwns2`: it catches
cases where the final alignment would still be identical (or trivially
differ) but the DP took a different intermediate path — useful because many
bugs only show up as intermediate-score drift that happens to cancel out.

Regenerate by temporarily adding an `fprintf` after the `pscore = ...`
assignment in the `treebase()` function in `mafft-upstream/core/disttbfast.c`
(around line 2108, immediately after `nlen[m1] = 0.5 * ...`):

```c
fprintf(stderr, "RDBG %d %d %d %d %.1f\n", l, clus1, clus2,
        (int)strlen(mseq1[0]), pscore);
```

Then rebuild C and run:

```bash
make -C mafft-upstream/core
MAFFT_BINARIES=$PWD/mafft-upstream/binaries \
  $PWD/mafft-upstream/scripts/mafft --nofft $PWD/mafft-upstream/test/sample \
  > /dev/null 2>&1 \
  | sed 's/.*RDBG/RDBG/g' \
  | LC_NUMERIC=C awk '/^RDBG/ {printf "%s %s %s %s %s %.1f\n", $1,$2,$3,$4,$5,$6}' \
  > crates/mafft-core/tests/fixtures/sample.nwns2.steps
```

Remember to revert the `fprintf` and rebuild C afterward, or the change will
leak into future regenerations of `sample.nwns2`.
