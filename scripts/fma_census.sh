#!/usr/bin/env bash
# Count fused multiply-add instructions in a C MAFFT build.
#
# C MAFFT 7.526 is not bit-reproducible across CPU architectures because the
# compiler decides whether `a * b + c` is contracted into a single-rounding
# FMA. This script makes the decision visible for a given install: it prints,
# for each of the three alignment engines (`disttbfast`, `dvtditr`, `tbfast`),
# the number of fused instructions in the disassembly. A non-zero count means
# the build fuses (mirror it with rust-MAFFT's default on aarch64, or
# `--features fp-contract-fma`); zero means it does not (mirror it with
# `--features fp-contract-none`, the default on non-aarch64 targets). See
# `crates/mafft-types/src/fp.rs`.
#
# Usage: scripts/fma_census.sh <dir-with-C-binaries>
#   e.g. scripts/fma_census.sh ~/.local/libexec/mafft          # manual install
#        scripts/fma_census.sh "$CONDA_PREFIX/libexec/mafft"   # bioconda
#        scripts/fma_census.sh mafft-upstream/binaries          # in-tree build
#
# macOS: uses `otool -tv` and counts AArch64 `fmadd/fmsub/fnmadd/fnmsub`
#        (scalar) and `fmla/fmls` (vector) mnemonics.
# Linux: uses `objdump -d` and counts x86-64 `vfmadd*/vfmsub*/vfnmadd*/
#        vfnmsub*` (FMA3, only present with -mfma / -march=native) and, on
#        aarch64 Linux, the same AArch64 mnemonics as macOS.
set -euo pipefail

if [ $# -ne 1 ] || [ ! -d "$1" ]; then
    echo "usage: $0 <dir containing disttbfast dvtditr tbfast>" >&2
    exit 2
fi
dir=$1

case "$(uname -s)" in
    Darwin)
        disasm() { otool -tv "$1"; }
        # Vector forms carry an arrangement suffix (`fmla.2d`, `fmla.d`).
        pattern='^[[:space:]]*[0-9a-fx]+[[:space:]]+(fmadd|fmsub|fnmadd|fnmsub|fmla|fmls)(\.[0-9a-z]+)?[[:space:]]'
        ;;
    Linux)
        disasm() { objdump -d "$1"; }
        case "$(uname -m)" in
            aarch64|arm64) pattern='[[:space:]](fmadd|fmsub|fnmadd|fnmsub|fmla|fmls)(\.[0-9a-z]+)?[[:space:]]' ;;
            *)             pattern='[[:space:]](vfmadd|vfmsub|vfnmadd|vfnmsub)[0-9]*[a-z]*[[:space:]]' ;;
        esac
        ;;
    *)
        echo "unsupported OS: $(uname -s)" >&2
        exit 2
        ;;
esac

printf '%-12s %10s %10s\n' binary total_insns fused_ops
for bin in disttbfast dvtditr tbfast; do
    path="$dir/$bin"
    if [ ! -x "$path" ]; then
        printf '%-12s %10s %10s\n' "$bin" missing missing
        continue
    fi
    listing=$(disasm "$path")
    total=$(printf '%s\n' "$listing" | grep -cE '^[[:space:]]*[0-9a-fx]+[[:space:]:]' || true)
    fused=$(printf '%s\n' "$listing" | grep -cE "$pattern" || true)
    printf '%-12s %10s %10s\n' "$bin" "$total" "$fused"
done
