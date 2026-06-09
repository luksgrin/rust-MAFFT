#!/usr/bin/env python3
"""BALIBASE byte-parity sweep for rust-MAFFT vs C MAFFT 7.526.

Runs both `mafft` (C) and `mafft-rs` on every unaligned `.tfa` input in a
BALIBASE distribution, diffs the FASTA outputs, and reports per-test byte
divergence counts plus an aggregate pass/fail summary.

The 36-seq sample in `mafft-upstream/test/sample` is too narrow to expose
all latent divergences (some show up only on certain length / divergence
profiles). BALIBASE 4 has ~218 curated test sets across reference families
RV11 / RV12 / RV20 / RV30 / RV40 / RV50, covering a wide range of
sequence counts and divergence levels — running mafft-rs against C on
every one of them is the highest-confidence way to confirm production
parity.

USAGE
-----

    scripts/balibase_parity.py BALIBASE_DIR [--modes MODE [MODE ...]] \\
        [--filter SUBSTR] [--limit N] [--timeout SECS] [--quiet] \\
        [--mafft PATH] [--mafft-rs PATH] [--output TSV]

Each MODE is a quoted string of mafft flags (e.g. "--maxiterate 100",
"--localpair", ""). The empty string means default mode (FFT-NS-2).

EXAMPLES
--------

    # Default mode (FFT-NS-2) across all tests:
    scripts/balibase_parity.py /path/to/BAliBASE/

    # Multiple modes:
    scripts/balibase_parity.py /path/to/BAliBASE/ \\
        --modes "" "--maxiterate 100" "--localpair" "--globalpair"

    # Smoke-test on a single reference set, with a hard cap:
    scripts/balibase_parity.py /path/to/BAliBASE/ \\
        --filter RV11 --limit 10 --timeout 60

    # Use a non-PATH C mafft binary:
    scripts/balibase_parity.py /path/to/BAliBASE/ \\
        --mafft /usr/local/bin/mafft

GETTING BALIBASE
----------------

BALIBASE 4 is distributed at https://www.lbgi.fr/balibase/. The structure
expected by this script is the unpacked archive, which contains directories
RV11 / RV12 / etc., each holding `.tfa` (unaligned input) and `.msf`
(reference) files. Only `.tfa` is used here — we don't compare against the
reference alignment, only against C MAFFT's output for the same input.

EXIT CODE
---------

Returns 0 if every (mode, test) pair has diff_lines == 0 or both
implementations errored identically. Returns 1 if any test has a real
parity gap. Useful for CI integration.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


# Default modes to sweep when --modes isn't supplied. The leading empty
# string is the default (no flags) FFT-NS-2 path.
DEFAULT_MODES = [""]


@dataclass
class Result:
    mode: str
    test_path: str
    nseq: int
    c_lines: int          # `wc -l` of C output; -1 = error
    r_lines: int          # `wc -l` of mafft-rs output; -1 = error
    diff_lines: int       # diff line count; -1 = either side errored
    status: str           # "match" | "diverge" | "c_error" | "rs_error" | "both_error" | "timeout"
    elapsed_c: float = 0.0
    elapsed_r: float = 0.0
    notes: str = ""


def count_fasta_seqs(path: Path) -> int:
    """Quick count of '>' header lines in a FASTA file."""
    try:
        with path.open("rb") as f:
            return sum(1 for line in f if line.startswith(b">"))
    except OSError:
        return 0


def run_mafft(binary: str, mode_args: list[str], input_path: Path,
              timeout: int) -> tuple[str | None, float, str]:
    """Run a MAFFT-style binary on `input_path` with `mode_args`.

    Returns (stdout_string_or_None, elapsed_seconds, error_summary).
    `stdout_string` is None if the run failed; in that case
    `error_summary` describes why ("timeout" / "exit_N" / "no_output").
    """
    argv = [binary, *mode_args, str(input_path)]
    start = time.perf_counter()
    try:
        proc = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return None, time.perf_counter() - start, "timeout"
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        return None, elapsed, f"exit_{proc.returncode}"
    if not proc.stdout:
        return None, elapsed, "no_output"
    return proc.stdout, elapsed, ""


def run_one(args, test: Path, mode: str) -> Result:
    """Run both implementations on a single test in a single mode."""
    mode_args = [] if mode == "" else mode.split()
    nseq = count_fasta_seqs(test)

    c_out, c_t, c_err = run_mafft(args.mafft, mode_args, test, args.timeout)
    r_out, r_t, r_err = run_mafft(args.mafft_rs, mode_args, test, args.timeout)

    if c_err and r_err:
        return Result(mode, str(test), nseq, -1, -1, -1, "both_error",
                      c_t, r_t, f"c:{c_err} rs:{r_err}")
    if c_err:
        return Result(mode, str(test), nseq, -1,
                      r_out.count("\n") if r_out else -1, -1,
                      "c_error", c_t, r_t, f"c:{c_err}")
    if r_err:
        return Result(mode, str(test), nseq, c_out.count("\n"), -1, -1,
                      "rs_error", c_t, r_t, f"rs:{r_err}")

    c_lines = c_out.count("\n")
    r_lines = r_out.count("\n")
    # Cheap-equality fast path: identical strings → diff is 0.
    if c_out == r_out:
        return Result(mode, str(test), nseq, c_lines, r_lines, 0, "match", c_t, r_t)

    # Compute diff line count via /tmp files (mirrors what users would do
    # interactively with `diff`).
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "c.fa").write_text(c_out)
        (td / "r.fa").write_text(r_out)
        diff = subprocess.run(
            ["diff", str(td / "c.fa"), str(td / "r.fa")],
            capture_output=True, text=True, check=False,
        )
    diff_lines = diff.stdout.count("\n")
    return Result(mode, str(test), nseq, c_lines, r_lines, diff_lines,
                  "diverge", c_t, r_t)


def emit_progress(r: Result, idx: int, total: int) -> None:
    status_icon = {
        "match": "✓",
        "diverge": "✗",
        "c_error": "C!",
        "rs_error": "R!",
        "both_error": "B!",
        "timeout": "T!",
    }.get(r.status, "?")
    name = r.test_path.split("/")[-1]
    mode = r.mode or "default"
    extra = f"diff={r.diff_lines}" if r.status == "diverge" else r.notes
    print(f"[{idx:4d}/{total}] {status_icon} {mode:<25} {name:<28} "
          f"nseq={r.nseq:<4} c={r.elapsed_c:5.2f}s r={r.elapsed_r:5.2f}s "
          f"{extra}", file=sys.stderr)


def write_tsv(rows: list[Result], path: Path) -> None:
    with path.open("w") as f:
        cols = ["mode", "test_path", "nseq", "c_lines", "r_lines",
                "diff_lines", "status", "elapsed_c", "elapsed_r", "notes"]
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write(f"{r.mode}\t{r.test_path}\t{r.nseq}\t{r.c_lines}\t"
                    f"{r.r_lines}\t{r.diff_lines}\t{r.status}\t"
                    f"{r.elapsed_c:.3f}\t{r.elapsed_r:.3f}\t{r.notes}\n")


def print_summary(rows: list[Result]) -> bool:
    """Print per-mode summary table. Return True if all tests matched."""
    by_mode: dict[str, list[Result]] = {}
    for r in rows:
        by_mode.setdefault(r.mode, []).append(r)
    all_clean = True
    print()
    print(f"{'Mode':<28} {'tests':>5} {'match':>6} {'diverge':>8} "
          f"{'c_err':>6} {'rs_err':>6} {'both':>5} {'time-out':>8}")
    print("-" * 80)
    for mode in sorted(by_mode):
        rs = by_mode[mode]
        counts = {"match": 0, "diverge": 0, "c_error": 0,
                  "rs_error": 0, "both_error": 0, "timeout": 0}
        for r in rs:
            counts[r.status] = counts.get(r.status, 0) + 1
        if counts["diverge"] or counts["rs_error"]:
            all_clean = False
        label = mode or "(default)"
        print(f"{label:<28} {len(rs):>5} {counts['match']:>6} "
              f"{counts['diverge']:>8} {counts['c_error']:>6} "
              f"{counts['rs_error']:>6} {counts['both_error']:>5} "
              f"{counts['timeout']:>8}")

    # Per-mode diverge / rs-error listing.
    diverges = [r for r in rows
                if r.status in ("diverge", "rs_error", "timeout")]
    if diverges:
        print()
        print("=== Non-matching tests (real parity gaps) ===")
        for r in diverges:
            tag = (f"diff={r.diff_lines}" if r.status == "diverge"
                   else r.status.upper())
            print(f"  {tag:<14} mode={r.mode or '(default)':<22} "
                  f"{r.test_path}  ({r.notes})")
    return all_clean


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="BALIBASE byte-parity sweep: mafft-rs vs C MAFFT",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("USAGE")[0],  # show only the intro
    )
    p.add_argument("balibase_dir", type=Path,
                   help="path to unpacked BALIBASE root (contains RV11/, "
                        "RV12/, etc. with .tfa inputs)")
    p.add_argument("--modes", nargs="*", default=DEFAULT_MODES,
                   help='mafft mode flag strings (e.g. "" "--maxiterate 100"). '
                        'Default: "" (FFT-NS-2 default mode).')
    p.add_argument("--filter", default="",
                   help="substring filter on input path (e.g. 'RV11', 'BB1')")
    p.add_argument("--pattern", default="*.tfa",
                   help="glob pattern for inputs under BALIBASE_DIR "
                        "(default '*.tfa'; for the Drive5 bench mirror use 'BB*')")
    p.add_argument("--limit", type=int, default=0,
                   help="run only first N tests after --filter (0 = all)")
    p.add_argument("--timeout", type=int, default=300,
                   help="per-alignment timeout in seconds (default 300)")
    p.add_argument("--mafft", default=os.environ.get("MAFFT") or shutil.which("mafft") or "mafft",
                   help="path to C mafft binary (default: $MAFFT or `which mafft`)")
    p.add_argument("--mafft-rs", default="./target/release/mafft-rs",
                   dest="mafft_rs",
                   help="path to mafft-rs binary (default: ./target/release/mafft-rs)")
    p.add_argument("--output", type=Path, default=Path("balibase_parity_results.tsv"),
                   help="TSV output path (default: balibase_parity_results.tsv)")
    p.add_argument("--quiet", action="store_true",
                   help="suppress per-test progress logging")
    return p.parse_args(argv)


def main() -> int:
    args = parse_args(sys.argv[1:])

    # Sanity checks
    if not args.balibase_dir.is_dir():
        print(f"error: BALIBASE dir not found: {args.balibase_dir}",
              file=sys.stderr)
        return 2
    if not Path(args.mafft_rs).is_file():
        print(f"error: mafft-rs binary not found: {args.mafft_rs}\n"
              "Build with: cargo build --release --bin mafft-rs",
              file=sys.stderr)
        return 2
    if not shutil.which(args.mafft) and not Path(args.mafft).is_file():
        print(f"error: C mafft not found: {args.mafft}", file=sys.stderr)
        return 2

    # Enumerate tests
    tests = sorted(p for p in args.balibase_dir.rglob(args.pattern) if p.is_file())
    if not tests:
        print(f"error: no files matching {args.pattern!r} under {args.balibase_dir}",
              file=sys.stderr)
        return 2
    if args.filter:
        tests = [t for t in tests if args.filter in str(t)]
    if args.limit > 0:
        tests = tests[:args.limit]

    total = len(tests) * len(args.modes)
    print(f"Running {total} alignments: {len(tests)} tests × "
          f"{len(args.modes)} mode(s)", file=sys.stderr)
    print(f"  mafft     : {args.mafft}", file=sys.stderr)
    print(f"  mafft-rs  : {args.mafft_rs}", file=sys.stderr)
    print(f"  timeout   : {args.timeout}s/test", file=sys.stderr)
    print(f"  output    : {args.output}", file=sys.stderr)
    print(file=sys.stderr)

    rows: list[Result] = []
    idx = 0
    for mode in args.modes:
        for test in tests:
            idx += 1
            row = run_one(args, test, mode)
            rows.append(row)
            if not args.quiet:
                emit_progress(row, idx, total)

    write_tsv(rows, args.output)
    print(f"\nWrote detailed results to {args.output}", file=sys.stderr)
    all_clean = print_summary(rows)
    return 0 if all_clean else 1


if __name__ == "__main__":
    sys.exit(main())
