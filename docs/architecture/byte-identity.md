# Byte-identity

rust-MAFFT's output matches C MAFFT 7.526 byte-for-byte across BAliBASE
3 (1930/1930 fixtures) for every supported mode, *when both are built on
the same platform*. This page explains why that's the bar, what it cost
to hit it, and why the "same platform" qualifier is not a hedge but a
property of the C reference itself.

## Why hold to byte-identity

MAFFT is the *de facto* MSA reference. Every downstream tool and
pipeline — phylogenetics, structure prediction, comparative genomics —
treats its output as ground truth. A Rust port that's "close" forces
every downstream pipeline to choose:

1. **Re-validate the rewrite** on their own data, accepting that
   alignment differences may cascade through downstream analyses.
2. **Maintain a parallel C dependency** so they can verify rust-MAFFT
   isn't introducing drift.

Both are non-starters for adoption. Byte-identity removes the choice:
the Rust binary is a drop-in replacement, no warm-up trust period
required.

## Parity status

### Modes covered

**Production-ready.** `mafft-rs` produces **byte-identical output to C MAFFT 7.526** across:

- All progressive modes: FFT-NS-2, NW-NS-2, FFT-NS-i, L-INS-1/i, G-INS-1/i, E-INS-1/i
- All scoring matrices: BLOSUM 30/45/50/62/80, JTT 100/200, TM 100/200 (FFT + NW + `--retree 1`)
- Tree variants: `--parttree`, `--dpparttree`, `--memsavetree`, `--youngestlinkage`, `--averagelinkage`, `--minimumlinkage`, `--mixedlinkage`, `--pileup`
- Iteration strategies: `--bestfirst`, `--simplehillclimbing`, `--skipiterate F` (both large- and small-F branches), `--oneiteration`
- Add-sequence modes: `--add`, `--addfragments`, `--keeplength` (including adversarial inputs with random insertions)
- Auxiliary outputs: `--distout`, `--scoreout`, `--nodeout`, `--mapout` / `--compactmapout`, `--treeout`
- DNA strand auto-detection: `--adjustdirection` (k-mer) and `--adjustdirectionaccurately` (DP)
- Refinement controls: `--allowshift`, `--unalignlevel`, `--maxambiguous`, `--minimumweight`, `--nwildcard`, `--leavegappyregion` / `--legacygappenalty`
- Memory strategy: `--memsave` / `--nomemsave` (Hirschberg DP), `--memsavetree`
- Seed alignment: `--seed FILE`, `--seedtable FILE` (single + multi-file)
- Fine-grained gap penalties: `--op`, `--ep`, `--exp`, `--shiftpenalty`, `--lop`, `--lep`, `--lexp`, `--gop`, `--gep`, `--gexp`
- Symlink dispatch via `argv[0]`: `linsi`, `ginsi`, `einsi`, `fftns`, `fftnsi`, `nwns`, `nwnsi`, `qinsi`, `xinsi`

Every progressive merge step matches in score and width; every refinement iteration converges to C's exact alignment.

### Parity matrix (36-seq protein sample, `mafft-upstream/test/sample`)

| Mode / flags                                | C width | Rust width | diff | Status |
|---------------------------------------------|---------|------------|------|--------|
| FFT-NS-2 (default)                          | 717     | 717        | 0    | ✓ byte-exact |
| NW-NS-2 (`--nofft`)                         | 717     | 717        | 0    | ✓ byte-exact |
| FFT-NS-i (`--maxiterate 100`)               | 721     | 721        | 0    | ✓ byte-exact |
| L-INS-1 (`--localpair --maxiterate 0`)      | 719     | 719        | 0    | ✓ byte-exact |
| L-INS-i (`linsi`)                           | 735     | 735        | 0    | ✓ byte-exact |
| G-INS-1 (`--globalpair --maxiterate 0`)     | 746     | 746        | 0    | ✓ byte-exact |
| G-INS-i (`ginsi`)                           | 737     | 737        | 0    | ✓ byte-exact |
| E-INS-1 (`--genafpair --maxiterate 0`)      | 729     | 729        | 0    | ✓ byte-exact |
| E-INS-i (`einsi`)                           | 729     | 729        | 0    | ✓ byte-exact |
| BLOSUM 30 / 45 / 50 / 62 / 80 (FFT + NW)    | match   | match      | 0    | ✓ byte-exact |
| JTT 100 / 200, TM 100 / 200 (FFT + NW)      | match   | match      | 0    | ✓ byte-exact |
| `--parttree`, `--dpparttree`, `--parttree --nofft`, `--parttree --reorder` | 752 | 752 | 0 | ✓ byte-exact |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6) | match | match | 0 | ✓ byte-exact |
| `--add` adversarial (16-existing + 3-added with random insertions) | match | match | 0 | ✓ byte-exact |
| `--pileup`, `--youngestlinkage`, `--averagelinkage`, `--minimumlinkage` | match | match | 0 | ✓ byte-exact |
| `--bestfirst` (FFT-NS-i / L-INS-i / G-INS-i)  | match | match      | 0    | ✓ byte-exact |
| `--skipiterate F` (F ∈ {0.05…1.0})         | match   | match      | 0    | ✓ byte-exact |
| `--oneiteration` (FFT-NS-2 + FFT-NS-i)      | match   | match      | 0    | ✓ byte-exact |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0    | ✓ byte-exact |
| `--allowshift --globalpair --maxiterate 1000` | 1060  | 1060       | 0    | ✓ byte-exact (all 386 BBaliBase 3) |
| `--adjustdirection`, `--adjustdirectionaccurately` (mixed forward/RC DNA) | match | match | 0 | ✓ byte-exact |
| `--memsavetree`, `--memsavetree --treeout`  | match   | match      | 0    | ✓ byte-exact |
| `--auto` (small/medium/large brackets)      | match   | match      | 0    | ✓ byte-exact |
| `--seed FILE`, `--seedtable FILE`           | match   | match      | 0    | ✓ byte-exact |
| `--treein`, `--treeout`, `--treein --treeout` | match | match      | 0    | ✓ byte-exact |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 0    | ✓ byte-exact |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 0‡   | ✓ byte-exact (needs `mxscarnamod`) |
| `--nuc`, `--amino` (forced sequence type)   | match   | match      | 0    | ✓ byte-exact |

‡ Q-INS-i requires the external `mxscarnamod` binary.

### DNA and input-handling parity (issue #1 integration)

| Input / flags | Status |
|---------------|--------|
| 120 × 1.4 kb CDS (`mtb_cds_120x1400.fa`), FFT-NS-i `--retree 2 --maxiterate 2` (1742 columns) | ✓ byte-exact |
| Same input, `--thread 1 --maxiterate 1/2/3/1000` (C's `athread` refinement rules) | ✓ byte-exact |
| First 8 sequences of the above, `--retree 2 --maxiterate 1000` (trailing `searchAnchors` region) | ✓ byte-exact |
| panaroo gene clusters (5 × 4 seqs, ~300–2300 bp), `--auto --adjustdirection --thread 1 --nuc` | ✓ byte-exact |
| DNA L-INS-1 / `--auto` on synthetic clusters (pair-phase gap scale ×3, `dndpre` offset 220) | ✓ 60/60 and 30/30 |
| Two-sequence input under `--maxiterate 1000` (C refines a pair exactly once, unweighted) | ✓ byte-exact |
| Input read path corpus (22 files × `default` / `--anysymbol` / `--localpair --maxiterate 0` / `--nuc` / `--add` modes, 71 cells incl. C's exit-1 cases) | ✓ 70/71 (see [What's still out of scope](index.md#whats-still-out-of-scope)) |

`--thread 1` and no `--thread` are *different* C code paths (`tditeration.c:1433`
selects `athread` on `nthread > 0`), with different branch order and
convergence rules; `mafft-rs` reproduces both byte-for-byte. C's own
`--thread N` for N ≥ 2 is nondeterministic, so byte-identity there is not
defined; `mafft-rs` stays deterministic at every thread count.

Residue case matches C MAFFT: nucleotide output is lowercased and protein
output uppercased at read time (C `io.c:1462-1467`, `io.c:1755`), driven by
the detected type or by `--nuc` / `--amino`; `--anysymbol` / `--preservecase`
keep the input's own case.

### BBaliBase 3 parity sweep

Sweep across BBaliBase 3 (RV11–RV50, 386 protein test sets) against MAFFT 7.526:

| Mode | Flags | Match | Total % |
|------|-------|-------|---------|
| FFT-NS-i | `--maxiterate 1000` | **386/386** | **100.0 %** |
| L-INS-i  | `--localpair --maxiterate 1000` | **386/386** | **100.0 %** |
| G-INS-i  | `--globalpair --maxiterate 1000` | **386/386** | **100.0 %** |
| E-INS-i  | `--genafpair --maxiterate 1000` | **386/386** | **100.0 %** |
| G-INS-i + allowshift | `--globalpair --maxiterate 1000 --allowshift` | **386/386** | **100.0 %** |

Total: **1930/1930 alignments byte-identical** to C MAFFT 7.526 across the BBaliBase 3 protein corpus.

## Byte-identity is per platform build

C MAFFT 7.526 is not bit-identical across CPU architectures. The
upstream Makefile compiles with `-O3` and no `-ffp-contract` setting,
which leaves the floating-point contraction policy to the compiler:

| C reference build | Compiler | Fused multiply-adds per binary |
|--|--|--|
| arm64 (Apple Silicon, `make` in-tree) | clang | on the order of a thousand — 1025 in the census that motivated this section; 882 scalar `fmadd`-family (1139 including vector `fmla`/`fmls`) with Apple clang 21, ~1140 in CI |
| x86-64 (`make` in-tree on Linux) | gcc | 0 |
| arm64 (Apple Silicon, `make CC=gcc-15` in-tree) | Homebrew GCC 15.2 | 0 |

clang fuses `a*b + c` into a single `fmadd` (one rounding); gcc on
x86-64 emits a multiply followed by an add (two roundings). Both are
valid C, and the resulting binaries disagree on some inputs. The policy
follows the compiler, not the architecture: GCC in an ISO mode such as
the Makefile's `-std=c99` does not contract even on arm64, so a GCC build
on Apple Silicon matches the `.nofma` fixtures (checked with
`MAFFT_FP_POLICY=nofma scripts/regen_policy_fixtures.sh … check`) and is
mirrored with `--features fp-contract-none`. The known
sentinel case is `--bl 50 --retree 2 --maxiterate 0` on the 36-sequence
`mafft-upstream/test/sample`:

| Platform | C MAFFT columns | rust-MAFFT columns (default policy) |
|--|--|--|
| arm64 / clang | 712 | 712 |
| x86-64 / gcc | 738 | 738 |

So "byte-identical to C" is always relative to a C build made on the
same platform, and rust-MAFFT tracks that platform rather than picking
one universal answer.

### The reference build

The reference is defined precisely: **the pinned `mafft-upstream` source,
built in-tree with the upstream Makefile's default flags
(`make -C mafft-upstream/core`), on the platform you run on.** No fixture
may be generated from a downloaded binary — a bioconda, Homebrew or
prebuilt-tarball `mafft` was compiled with someone else's toolchain and
flags, and it is not what the tests, the FFI cross-validation or CI compare
against.

### The `fp` policy module

`mafft_types::fp` centralises the decision:

- `mafft_types::fp::CONTRACTS_FMA` — a `const bool`. Defaults to `true`
  on `aarch64` (mirroring clang's behaviour there) and `false`
  elsewhere.
- `mafft_types::fp::fmadd(a, b, c)` — computes `a*b + c` with a single
  rounding when `CONTRACTS_FMA` is set, and as two separate operations
  otherwise. Every hot loop that has to match a contracted C expression
  goes through it rather than calling `f64::mul_add` directly.
- Cargo features `fp-contract-fma` / `fp-contract-none` override the
  default. They are defined on `mafft-types` and forwarded by every
  crate in the workspace, so `cargo test -p mafft-core --features
  fp-contract-none` builds the whole dependency graph under the x86-64
  policy on any host. (`--features` cannot be passed at the
  virtual-workspace root; use `-p <crate>`.)

### Fixture convention

Where a C-canonical fixture differs between the two policies it is
committed twice, as `<name>.fma` and `<name>.nofma`, and the test
helper `fixture_path_fp(name)` picks the one matching `CONTRACTS_FMA`.
Fixtures without a suffix are identical under both policies.

Each pair is described by a line in
`crates/mafft-core/tests/fixtures/policy_fixtures.tsv`
(`fixture_basename<TAB>input_relpath<TAB>mafft_args`) and regenerated by
`scripts/regen_policy_fixtures.sh`, which detects the policy from
`uname -m` (arm64/aarch64 => `.fma`, else `.nofma`; override with
`MAFFT_FP_POLICY=fma|nofma`):

```bash
make -C mafft-upstream/core
scripts/regen_policy_fixtures.sh mafft-upstream/scripts/mafft check   # diff against the checked-in variant
scripts/regen_policy_fixtures.sh mafft-upstream/scripts/mafft write   # overwrite this host's variant
```

A `.fma` fixture therefore comes from the in-tree build on an arm64 host
and a `.nofma` one from the in-tree build on an x86-64 host; never from
the other, and never from a binary you did not build from the pinned
source.

### What CI checks

- The `Build (x86-64)` and `Build (arm64)` jobs each build the C
  submodule in-tree and run the full suite (lib, integration and FFI
  `cross_validate_*` tests) against it, so the claim is checked on both
  reference platforms. Each job prints an FMA census of the C binaries
  it just built and runs the `--bl 50` sentinel above, diffing
  `mafft-rs` against the same-platform C `mafft` wrapper. It then runs
  `scripts/regen_policy_fixtures.sh ... check`, regenerating every
  policy-sensitive fixture from that in-tree C build and failing on any
  byte difference from the checked-in variant for the runner's policy
  (`.nofma` on x86-64, `.fma` on arm64) — so the committed fixtures are
  continuously re-derived from the reference definition rather than
  trusted.
- A `policy-cross-check` job runs the fixture-based (non-FFI) test
  binaries on the arm64 runner with `--features fp-contract-none`, which
  checks the `.nofma` fixtures from an arm64 host. The FFI tests are
  deliberately excluded there: they compare against the C compiled on
  the same runner, which contracts, so under `fp-contract-none` the two
  sides would be running different policies.

## Design decisions that fall out

### No global state

The C code uses ~400 `extern` globals. The Rust modules use owned
`ScoringContext`, `Topology`, `Profile` structs passed explicitly. The
engine is reentrant, embeddable, and trivially thread-safe at the
call-site boundary.

This wasn't done for ergonomic reasons — it was forced by byte-identity.
C MAFFT's globals get assigned in subtle orderings (e.g., `gapfaclocal`
depends on whether `--allowshift` is parsed before or after
`unalignlevel`). Replicating that in Rust required making the
"set-it-now" semantics explicit.

### Hand-ported Cooley-Tukey FFT

`mafft-fft/src/fft_c_compat.rs` is a hand port of MAFFT's `core/fft.c`,
matching the C rounding exactly. Off-the-shelf FFT libraries (rustfft,
realfft, FFTW) vary in butterfly grouping order, producing 1-ULP
correlation differences. On flat-landscape similarity matrices those
ULPs flip FFT-anchor selection — which propagates into different
alignment segmentation and ultimately different output bytes.

### FP order matters

In several hot loops the C compiler's contraction policy decides the
result: with clang on arm64, `a*b + c` becomes a single fused
multiply-add (one rounding); with gcc on x86-64 it stays two operations
(two roundings). rust-MAFFT reproduces whichever the host's C build
does through `mafft_types::fp::fmadd`, which fuses only when
`CONTRACTS_FMA` is set (see [above](#byte-identity-is-per-platform-build)).
One concrete case (`calcW` in `mafft-tree/src/weighting.rs`), matching
clang's `fmadd d2, b, c, d3; fmadd d0, a, b, d2` sequence:

```rust
// a*b + (b*c + a*c), fused exactly where clang fuses it; under the
// x86-64 policy `fmadd` degrades to the unfused mul + add gcc emits.
let s = fmadd(a, b, fmadd(b, c, a * c));
```

Doing this as a naive `a*b + b*c + a*c` on arm64 produces visibly
different weights on `BB40043`, `BB30018`, `BB40010`, `BB30010`,
`BB40004` — and doing it fused on x86-64 would drift from the gcc build
the same way.

### Tie-break ordering

Tied scores in Needleman-Wunsch / Smith-Waterman traceback are resolved
in C via `if (new >= best)` (note the `>=` — equality favors the new
candidate). Several constraint-mode scanners use a different inequality
(`if (new > best)`). Mixing them up flips one trace direction on one
column, which propagates downstream. We mirror each C usage site
exactly.

### Per-step weight normalization

C MAFFT's `Falign` (the FFT-anchored DP entry) expects the per-group
weight vectors to sum to 1.0. The Rust progressive-alignment caller
normalizes before invocation. Skipping this on a single merge step
produces an off-by-one column on `BB30013`.

## Cross-validation harness

The byte-identity bar is held in place by 82 FFI test functions across
15 files (`crates/*/tests/cross_validate*.rs`; 99 including the three
forensic binaries) that compile MAFFT's C source in-tree via
`mafft-c-bindings` and call both implementations on the same inputs,
comparing outputs byte-for-byte. Because the C is
compiled on the same machine, these tests check Rust's default policy
against the matching C build on whichever platform they run.

This is documented in detail on the
[Cross-validation page](cross-validation.md).

## What this means for changes

Any change to `mafft-core` / `mafft-align` / `mafft-fft` / `mafft-tree`
/ `mafft-scoring` MUST keep the test suite green. New behaviour goes
behind a feature flag or a new function — not as a tweak to an existing
hot path. If you're tempted to "clean up" a slightly weird-looking
arithmetic in the engine, check the cross-validate test for that
function first; it's probably weird because the C output requires it.

## What this DOESN'T mean

- Algorithmic improvements are still possible — they just need their
  own modes / flags, and the byte-identical existing modes must keep
  working.
- Performance optimisation is fine wherever it doesn't change the
  observable output. We've shaved cycles via pre-computed boundary
  tables and LTO+fat codegen — none of which alter the alignment.
  Collapsing arithmetic into `f64::mul_add` is *not* in this category:
  it changes rounding, so it must go through `fp::fmadd` and be
  justified by the C assembly on the matching platform.
- The Python and CLI surfaces can evolve independently of the engine.
  Pretty-printing, progress callbacks, alternate output formats — all
  on the table.
