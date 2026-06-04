#!/usr/bin/env bash
# Render the `mafft-rs` CLI reference as a markdown file by running
# `mafft-rs --help` and wrapping the output in a fenced block.
#
# Usage:
#   ./scripts/gen-cli-reference.sh > docs/cli-reference.md
#
# The CI workflow (`.github/workflows/docs.yml`) invokes this after
# `cargo build --release -p mafft-rs`, before `mkdocs build`. Locally,
# you can do the same:
#
#   cargo build --release -p mafft-rs
#   ./scripts/gen-cli-reference.sh > docs/cli-reference.md
#   mkdocs serve

set -euo pipefail

# Resolve the binary. Honor MAFFT_RS_BIN if set (CI uses this), else
# fall back to the workspace's release build.
if [[ -n "${MAFFT_RS_BIN:-}" ]]; then
    BIN="${MAFFT_RS_BIN}"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
    REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
    BIN="${REPO_ROOT}/target/release/mafft-rs"
fi

if [[ ! -x "${BIN}" ]]; then
    echo "error: ${BIN} not found or not executable" >&2
    echo "       run: cargo build --release -p mafft-rs" >&2
    exit 1
fi

VERSION="$("${BIN}" --version | head -1)"

cat <<EOF
# CLI Reference

Auto-generated from \`mafft-rs --help\` on every docs build. The
authoritative source for every flag is the CLI itself — run
\`mafft-rs --help\` locally for the version you have installed.

**Generated for**: \`${VERSION}\`

## Usage

\`\`\`text
EOF

"${BIN}" --help

cat <<'EOF'
```

## Notes on flag semantics

- `--maxiterate 0` disables iterative refinement (FFT-NS-2 default).
- `--maxiterate N` (N > 0) enables refinement up to N cycles per
  strategy. INS-i strategies default to 1000 if `--maxiterate` is
  unset.
- `--quiet` suppresses the per-step progress noise to stderr; the
  alignment itself still writes to stdout.
- `--nofft` forces NW/SW pairwise alignment for every merge step —
  exact behavior, but slower than FFT-anchored DP.
- `--pdbidlist` and `--pdbfilelist` print a "temporarily unavailable,
  2018/Dec." message and exit 0, matching the C MAFFT 7.526 behavior
  verbatim (these flags were disabled upstream).

## Adding new flags

CLI flags live in `crates/mafft-bin/src/lib.rs` as fields on the
`Args` struct (parsed via `clap` derive macros). After adding a flag,
re-run this script to regenerate the reference; the next docs build
will pick it up automatically.
EOF
