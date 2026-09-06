#!/usr/bin/env bash
# Regenerate or check the policy-sensitive test fixtures against the in-tree
# C MAFFT build.
#
# C MAFFT 7.526 is not bit-reproducible across CPU architectures: the
# compiler decides whether `a * b + c` is contracted into a single-rounding
# fused multiply-add (clang on arm64 does, gcc on baseline x86-64 does not),
# and on some inputs that changes the alignment. Fixtures whose bytes depend
# on it are committed twice, as `<name>.fma` and `<name>.nofma`
# (see crates/mafft-core/tests/fixtures/README.md).
#
# The reference for BOTH variants is the same thing: the pinned
# `mafft-upstream` source, built in-tree with the upstream Makefile's default
# flags (`make -C mafft-upstream/core`), on the platform you run on. Which
# variant that build corresponds to is decided by the host architecture
# (arm64/aarch64 => .fma, anything else => .nofma). No fixture may be
# generated from a downloaded binary.
#
# Usage:
#   scripts/regen_policy_fixtures.sh <path/to/mafft-upstream/scripts/mafft> check
#   scripts/regen_policy_fixtures.sh <path/to/mafft-upstream/scripts/mafft> write
#
#   check  regenerate every manifest entry and diff (unified) against the
#          checked-in variant for this host's policy; exit 1 on any
#          difference or missing fixture.
#   write  overwrite the variant for this host's policy with the freshly
#          generated output (the other variant is left untouched).
#
# Environment:
#   MAFFT_FP_POLICY=fma|nofma  override the policy detected from `uname -m`.
#   MAFFT_BINARIES             where the wrapper finds its engines; defaults
#                              to `<wrapper dir>/../binaries` (populated by
#                              `make -C mafft-upstream/core`).
#
# Manifests: every `crates/*/tests/fixtures/policy_fixtures.tsv`, with
# tab-separated columns `fixture_basename  input_relpath  mafft_args`
# (`#` comments and blank lines ignored; input_relpath is relative to the
# repository root; mafft_args are whitespace-split, no quoting). The script
# also fails if a `.fma`/`.nofma` file exists under `crates/*/tests/fixtures`
# that no manifest describes, so a new policy-sensitive fixture cannot be
# added without saying how it is produced.
#
# Portable to bash 3.2 (macOS) and GNU bash (Linux).
set -euo pipefail

usage() {
    echo "usage: $0 <path/to/mafft-upstream/scripts/mafft> <check|write>" >&2
    exit 2
}

[ $# -eq 2 ] || usage
wrapper=$1
mode=$2
case "$mode" in
    check|write) ;;
    *) echo "error: mode must be 'check' or 'write', got '$mode'" >&2; usage ;;
esac
if [ ! -x "$wrapper" ]; then
    echo "error: C mafft wrapper '$wrapper' is not executable (run 'make -C mafft-upstream/core' first)" >&2
    exit 2
fi

script_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)
wrapper=$(cd "$(dirname "$wrapper")" && pwd)/$(basename "$wrapper")

# The in-tree wrapper reads MAFFT_BINARIES; the Makefile's `all` target
# copies the engines into ../binaries relative to scripts/.
if [ -z "${MAFFT_BINARIES:-}" ]; then
    candidate=$(cd "$(dirname "$wrapper")/.." && pwd)/binaries
    if [ -x "$candidate/disttbfast" ]; then
        export MAFFT_BINARIES=$candidate
    else
        echo "error: MAFFT_BINARIES is unset and '$candidate/disttbfast' does not exist" >&2
        exit 2
    fi
fi

# Policy: explicit override, else host architecture.
arch=$(uname -m)
if [ -n "${MAFFT_FP_POLICY:-}" ]; then
    policy=$MAFFT_FP_POLICY
    policy_source="MAFFT_FP_POLICY"
else
    case "$arch" in
        arm64|aarch64) policy=fma ;;
        *)             policy=nofma ;;
    esac
    policy_source="uname -m = $arch"
fi
case "$policy" in
    fma|nofma) ;;
    *) echo "error: MAFFT_FP_POLICY must be 'fma' or 'nofma', got '$policy'" >&2; exit 2 ;;
esac
other=nofma
[ "$policy" = nofma ] && other=fma

echo "== policy-sensitive fixtures: mode=$mode policy=.$policy ($policy_source)"
echo "   C wrapper:      $wrapper"
echo "   MAFFT_BINARIES: $MAFFT_BINARIES"

tmp=$(mktemp -d "${TMPDIR:-/tmp}/regen_policy_fixtures.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

manifests=$(find "$repo_root/crates" -path '*/tests/fixtures/policy_fixtures.tsv' -type f | LC_ALL=C sort)
if [ -z "$manifests" ]; then
    echo "error: no crates/*/tests/fixtures/policy_fixtures.tsv manifest found" >&2
    exit 2
fi

status=0
listed=$tmp/listed
: > "$listed"

for manifest in $manifests; do
    fixdir=$(dirname "$manifest")
    # `read` with IFS=tab keeps the args column intact; a trailing field
    # without a newline is still delivered thanks to the `|| [ -n ... ]`.
    while IFS=$'\t' read -r name input args || [ -n "${name:-}" ]; do
        case "${name:-}" in ''|'#'*) continue ;; esac
        if [ -z "${input:-}" ]; then
            echo "error: $manifest: entry '$name' has no input column" >&2
            status=1; continue
        fi
        expected=$fixdir/$name.$policy
        echo "$expected" >> "$listed"
        echo "$fixdir/$name.$other" >> "$listed"
        input_abs=$repo_root/$input
        if [ ! -f "$input_abs" ]; then
            echo "error: $name: input '$input' not found under $repo_root" >&2
            status=1; continue
        fi
        actual=$tmp/$name.$policy
        echo "-- $name.$policy  <=  mafft ${args:-} $input"
        # shellcheck disable=SC2086  # args are meant to be word-split
        if ! "$wrapper" ${args:-} "$input_abs" > "$actual"; then
            echo "error: $name: C mafft failed" >&2
            status=1; continue
        fi
        case "$mode" in
            write)
                cp "$actual" "$expected"
                echo "   wrote $(wc -c < "$expected" | tr -d ' ') bytes to ${expected#"$repo_root"/}"
                ;;
            check)
                if [ ! -f "$expected" ]; then
                    echo "error: $name: missing checked-in fixture ${expected#"$repo_root"/}" >&2
                    status=1; continue
                fi
                if cmp -s "$expected" "$actual"; then
                    echo "   OK: identical to ${expected#"$repo_root"/}"
                else
                    echo "error: $name: in-tree C output differs from ${expected#"$repo_root"/}" >&2
                    diff -u -L "${expected#"$repo_root"/} (checked in)" \
                            -L "$name.$policy (in-tree C, $arch)" \
                            "$expected" "$actual" || true
                    if [ -f "$fixdir/$name.$other" ] && cmp -s "$fixdir/$name.$other" "$actual"; then
                        echo "note: the output matches $name.$other instead; is the policy detection wrong for this host?" >&2
                    fi
                    status=1
                fi
                ;;
        esac
    done < "$manifest"
done

# Every committed variant must be described by a manifest entry.
unlisted=0
for f in $(find "$repo_root/crates" -path '*/tests/fixtures/*' -type f \( -name '*.fma' -o -name '*.nofma' \) | LC_ALL=C sort); do
    if ! grep -qxF "$f" "$listed"; then
        echo "error: ${f#"$repo_root"/} is not described by any policy_fixtures.tsv manifest" >&2
        unlisted=1
    fi
done
[ $unlisted -eq 0 ] || status=1

if [ $status -eq 0 ]; then
    echo "== policy-sensitive fixtures: $mode OK (.$policy)"
else
    echo "== policy-sensitive fixtures: $mode FAILED (.$policy)" >&2
fi
exit $status
