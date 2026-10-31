#!/usr/bin/env python3
"""Run inputs through a sanitizer-instrumented C MAFFT and through mafft-rs (issue #27).

Each run is executed twice: once with the C reference (the pinned
`mafft-upstream` compiled with `-fsanitize=address,undefined`, driven by
its own `mafft` wrapper script) and once with the release `mafft-rs`
binary. Per run the harness records:

* C: exit status, wall time, every AddressSanitizer / UndefinedBehaviorSanitizer
  report (counted from the sanitizer log files, so reports from engines whose
  stderr the wrapper discards are not lost), the first report line with its
  first `mafft-upstream` stack frame, signals reported by the shell
  (`Segmentation fault`, `Abort trap`, ...) and a short stderr excerpt;
* Rust: exit status, wall time, outcome (`ok`, `ok-empty`, `error`, `panic`,
  `signal`, `timeout`) and a short stderr excerpt;
* `identical`: whether stdout is byte-identical when both exit 0.

Rows are appended to `benchmarks/results/memsafety_<date>.tsv` as runs
finish; re-running skips every (corpus, case, flags) already present, so an
interrupted sweep resumes where it stopped. Full stderr and sanitizer logs
go to `target/memsafety/logs/`. The harness renices itself to 19.

PREREQUISITES
-------------

    # sanitized C build, installed into target/asan/mafft-bin/
    nice -n 19 make -C mafft-upstream/core -j6 CC=clang STRIP=true \\
        CFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' \\
        LDFLAGS='-fsanitize=address,undefined' \\
        PREFIX=$PWD/target/asan/mafft-bin install
    make -C mafft-upstream/core clean
    cargo build --release -p mafft-rs
    benchmarks/memsafety_gen.py          # generated sets under target/memsafety/

USAGE
-----

    benchmarks/memsafety.py                                  # all small corpora
    benchmarks/memsafety.py --corpus input_handling sample
    benchmarks/memsafety.py --corpus balibase --balibase DIR # the deferred sweep
    benchmarks/memsafety.py --summary                        # summarise the TSV

Corpora: input_handling, sample, clusters, adversarial, sample_add,
mutants, cases (`benchmarks/memsafety_cases/`), balibase (not in the
default set).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FIX = REPO / "crates" / "mafft-bin" / "tests" / "fixtures"
GEN = REPO / "target" / "memsafety"
LOGS = GEN / "logs"
CASES = REPO / "benchmarks" / "memsafety_cases"
ASAN_PREFIX = REPO / "target" / "asan" / "mafft-bin"

COLUMNS = [
    "corpus", "case", "flags", "input",
    "c_exit", "c_secs", "c_asan", "c_ubsan", "c_signal", "c_diag", "c_stderr",
    "rs_exit", "rs_secs", "rs_outcome", "rs_stderr", "identical",
]
DEFAULT_CORPORA = ["input_handling", "sample", "clusters", "adversarial", "sample_add", "mutants", "cases"]

ASAN_RE = re.compile(r"ERROR: AddressSanitizer: .*")
UBSAN_RE = re.compile(r"\S+:\d+:\d+: runtime error: .*")
FRAME_RE = re.compile(r"#\d+ 0x[0-9a-f]+ in (\S+) (?:\S*/)?(\w+\.c:\d+)")
ACCESS_RE = re.compile(r"(READ|WRITE) of size \d+")
SIGNAL_RE = re.compile(r"(Segmentation fault|Bus error|Abort trap|Illegal instruction|Floating point exception|Killed)[^\n]*")

INPUT_HANDLING_MODES = {
    "default": "",
    "anysymbol": "--anysymbol",
    "linsi1": "--localpair --maxiterate 0",
    "nuc": "--nuc",
    "add": "--add {fix}/input_handling/add_new_unusual.fa",
}
SAMPLE_MODES = [
    "", "--retree 1", "--maxiterate 1000", "--globalpair --maxiterate 16",
    "--localpair --maxiterate 16", "--genafpair --maxiterate 16", "--parttree",
    "--dpparttree", "--memsave", "--adjustdirection", "--adjustdirectionaccurately",
    "--auto", "--thread 2", "--treeout", "--anysymbol",
]
SAMPLERNA_MODES = [
    "", "--nuc", "--globalpair --maxiterate 16", "--localpair --maxiterate 16",
    "--parttree", "--memsave", "--adjustdirection --nuc",
]
CLUSTER_MODES = ["--auto --adjustdirection --thread 1 --nuc", "--retree 2 --maxiterate 1000 --nuc"]
BALIBASE_MODES = ["", "--maxiterate 1000", "--localpair --maxiterate 1000"]


def rel(p: Path) -> str:
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def manifest_jobs(corpus: str, d: Path) -> list[tuple[str, str, str, Path]]:
    m = d / "manifest.tsv"
    if not m.exists():
        subprocess.run([sys.executable, str(REPO / "benchmarks" / "memsafety_gen.py")], check=True)
    jobs = []
    with open(m) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            jobs.append((corpus, row["case"], row["flags"].replace("{dir}", str(d)), d / row["input"]))
    return jobs


def build_jobs(corpora: list[str], balibase: Path | None) -> list[tuple[str, str, str, Path]]:
    jobs: list[tuple[str, str, str, Path]] = []
    for corpus in corpora:
        if corpus == "input_handling":
            for fa in sorted((FIX / "input_handling").glob("*.fa")):
                for mode, flags in INPUT_HANDLING_MODES.items():
                    jobs.append((corpus, fa.stem, flags.replace("{fix}", str(FIX)), fa))
        elif corpus == "sample":
            # Copied so `--treeout` never writes into the submodule.
            d = GEN / "upstream"
            d.mkdir(parents=True, exist_ok=True)
            for name, modes in (("sample", SAMPLE_MODES), ("samplerna", SAMPLERNA_MODES)):
                src = REPO / "mafft-upstream" / "test" / name
                dst = d / name
                shutil.copyfile(src, dst)
                jobs += [(corpus, name, flags, dst) for flags in modes]
        elif corpus == "clusters":
            fas = sorted((FIX / "panaroo_tiny_dna_clusters").glob("*.fa"))
            fas += sorted((FIX / "r2_dna_clusters").glob("*.fa"))
            fas.append(FIX / "mtb_cds_120x1400.fa")
            jobs += [(corpus, fa.stem, flags, fa) for fa in fas for flags in CLUSTER_MODES]
        elif corpus in ("adversarial", "sample_add", "mutants"):
            jobs += manifest_jobs(corpus, GEN / corpus)
        elif corpus == "cases":
            if (CASES / "manifest.tsv").exists():
                jobs += manifest_jobs(corpus, CASES)
        elif corpus == "balibase":
            if balibase is None:
                sys.exit("--corpus balibase needs --balibase DIR (the directory holding the RV* folders)")
            for fa in sorted(balibase.glob("RV*/*.tfa")) or sorted(balibase.glob("**/*.tfa")):
                jobs += [(corpus, fa.stem, flags, fa) for flags in BALIBASE_MODES]
        else:
            sys.exit(f"unknown corpus {corpus}")
    return jobs


def excerpt(text: str, n: int = 3, width: int = 240) -> str:
    """Last `n` meaningful stderr lines, without progress counters or sanitizer reports."""
    keep = []
    for line in re.split(r"[\r\n]+", text):
        s = line.strip()
        if not s or SIGNAL_RE.search(s) or "Sanitizer" in s or "runtime error:" in s:
            continue
        if re.match(r"^(STEP|\d+ / \d+|done\.|Making|Constructing|Progressive|Loading|Reallocating|nthread|stacksize|generating|tbfast|disttbfast)", s):
            continue
        keep.append(s)
    return " | ".join(keep[-n:])[:width].replace("\t", " ")


def run(cmd: list[str], env: dict[str, str], timeout: float, cwd: Path) -> tuple[int | None, float, bytes, str]:
    t0 = time.monotonic()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL,
                         env=env, cwd=cwd, start_new_session=True)
    try:
        out, err = p.communicate(timeout=timeout)
        code: int | None = p.returncode
    except subprocess.TimeoutExpired:
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        out, err = p.communicate()
        code = None
    return code, time.monotonic() - t0, out, err.decode("utf-8", "replace")


def run_job(job, args, logdir: Path) -> dict[str, str]:
    corpus, case, flags, inp = job
    key = hashlib.sha1(f"{corpus}\0{case}\0{flags}".encode()).hexdigest()[:12]
    stem = re.sub(r"[^\w.-]", "_", case)[:60] + "." + key
    logs = logdir / corpus
    logs.mkdir(parents=True, exist_ok=True)
    work = GEN / "work" / key
    work.mkdir(parents=True, exist_ok=True)
    argv = shlex.split(flags) + [str(inp)]

    san_prefix = logs / f"{stem}.san"
    for old in logs.glob(f"{stem}.san.*"):
        old.unlink()
    env = dict(os.environ)
    env["MAFFT_BINARIES"] = str(args.c_prefix / "libexec" / "mafft")
    env["ASAN_OPTIONS"] = f"detect_leaks=0:abort_on_error=0:log_path={san_prefix}"
    env["UBSAN_OPTIONS"] = f"print_stacktrace=1:log_path={san_prefix}"
    env["TMPDIR"] = str(work)
    c_code, c_secs, c_out, c_err = run([str(args.c_prefix / "bin" / "mafft")] + argv, env, args.timeout, work)

    san_text = "".join(f.read_text("utf-8", "replace") for f in sorted(logs.glob(f"{stem}.san.*")))
    san_text += c_err
    asan = ASAN_RE.findall(san_text)
    ubsan = sorted(set(UBSAN_RE.findall(san_text)))
    diag = ""
    if asan or ubsan:
        first = asan[0] if asan else ubsan[0]
        frame = FRAME_RE.search(san_text)
        diag = re.sub(r"\s+", " ", first)
        diag = re.sub(r" on address 0x[0-9a-f]+ at pc 0x[0-9a-f]+ bp 0x[0-9a-f]+ sp 0x[0-9a-f]+", "", diag)
        diag = diag.replace(str(REPO) + "/", "")
        access = ACCESS_RE.search(san_text) if asan else None
        if access:
            diag += f" ({access.group(0)})"
        if frame:
            diag += f" @ {frame.group(1)} {frame.group(2)}"
    sig = SIGNAL_RE.search(c_err)
    (logs / f"{stem}.c.err").write_text(c_err[-1_000_000:])

    renv = {k: v for k, v in os.environ.items() if not k.startswith(("ASAN", "UBSAN", "MAFFT_BINARIES"))}
    r_code, r_secs, r_out, r_err = run([str(args.rs)] + argv, renv, args.timeout, work)
    (logs / f"{stem}.rs.err").write_text(r_err[-1_000_000:])
    if r_code is None:
        outcome = "timeout"
    elif r_code < 0:
        outcome = f"signal{-r_code}"
    elif "panicked at" in r_err:
        outcome = "panic"
    elif r_code == 0:
        outcome = "ok" if r_out.strip() else "ok-empty"
    else:
        outcome = "error"
    identical = ""
    if c_code == 0 and r_code == 0:
        identical = "yes" if c_out == r_out else "no"
        if c_out != r_out:
            (logs / f"{stem}.c.out").write_bytes(c_out)
            (logs / f"{stem}.rs.out").write_bytes(r_out)
    shutil.rmtree(work, ignore_errors=True)

    return {
        "corpus": corpus, "case": case, "flags": flags.replace(str(REPO) + "/", ""), "input": rel(inp),
        "c_exit": "timeout" if c_code is None else str(c_code), "c_secs": f"{c_secs:.2f}",
        "c_asan": str(len(asan)), "c_ubsan": str(len(ubsan)),
        "c_signal": sig.group(1) if sig else "", "c_diag": diag[:300], "c_stderr": excerpt(c_err),
        "rs_exit": "timeout" if r_code is None else str(r_code), "rs_secs": f"{r_secs:.2f}",
        "rs_outcome": outcome, "rs_stderr": excerpt(r_err), "identical": identical,
    }


def summarise(path: Path) -> None:
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    by: dict[str, dict[str, int]] = {}
    for r in rows:
        s = by.setdefault(r["corpus"], dict(runs=0, c_ok=0, c_fail=0, san=0, sig=0, rs_ok=0, rs_err=0,
                                             rs_bad=0, same=0, diff=0, agree=0))
        s["runs"] += 1
        s["c_ok" if r["c_exit"] == "0" else "c_fail"] += 1
        s["san"] += int(r["c_asan"] != "0" or r["c_ubsan"] != "0")
        s["sig"] += int(bool(r["c_signal"]))
        s["rs_ok" if r["rs_outcome"].startswith("ok") else "rs_err" if r["rs_outcome"] == "error" else "rs_bad"] += 1
        s["same"] += int(r["identical"] == "yes")
        s["diff"] += int(r["identical"] == "no")
        s["agree"] += int((r["c_exit"] == "0") == (r["rs_exit"] == "0"))
    cols = list(next(iter(by.values())).keys()) if by else []
    print("corpus\t" + "\t".join(cols))
    for c, s in by.items():
        print(c + "\t" + "\t".join(str(s[k]) for k in cols))
    print("\nruns with a sanitizer report or a signal:")
    for r in rows:
        if r["c_asan"] != "0" or r["c_ubsan"] != "0" or r["c_signal"]:
            print(f"  {r['corpus']}/{r['case']} [{r['flags']}] C={r['c_exit']} {r['c_signal']} {r['c_diag']}"
                  f" || rs={r['rs_outcome']}({r['rs_exit']}) {r['rs_stderr'][:120]}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", nargs="+", default=DEFAULT_CORPORA)
    ap.add_argument("--out", type=Path, default=REPO / "benchmarks" / "results" / f"memsafety_{dt.date.today()}.tsv")
    ap.add_argument("--timeout", type=float, default=120.0)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--filter", default="", help="regex on the case name")
    ap.add_argument("--c-prefix", type=Path, default=ASAN_PREFIX, help="install prefix of the sanitized C build")
    ap.add_argument("--rs", type=Path, default=REPO / "target" / "release" / "mafft-rs")
    ap.add_argument("--balibase", type=Path, default=None)
    ap.add_argument("--summary", action="store_true", help="only summarise --out")
    args = ap.parse_args()

    if args.summary:
        summarise(args.out)
        return
    os.nice(19)
    for p in (args.c_prefix / "bin" / "mafft", args.rs):
        if not p.exists():
            sys.exit(f"missing {p} (see PREREQUISITES in --help)")

    done: set[tuple[str, str, str]] = set()
    if args.out.exists():
        with open(args.out) as f:
            done = {(r["corpus"], r["case"], r["flags"]) for r in csv.DictReader(f, delimiter="\t")}
    jobs = [j for j in build_jobs(args.corpus, args.balibase)
            if (j[0], j[1], j[2].replace(str(REPO) + "/", "")) not in done and re.search(args.filter, j[1])]
    print(f"{len(jobs)} runs to do ({len(done)} already in {rel(args.out)})", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    new = not args.out.exists()
    lock = threading.Lock()
    with open(args.out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
            f.flush()
        count = [0]

        def one(job):
            row = run_job(job, args, LOGS)
            with lock:
                w.writerow(row)
                f.flush()
                count[0] += 1
                flag = " SANITIZER" if row["c_diag"] else ""
                print(f"[{count[0]}/{len(jobs)}] {row['corpus']}/{row['case']} [{row['flags']}] "
                      f"C={row['c_exit']} rs={row['rs_outcome']} same={row['identical'] or '-'}{flag}",
                      file=sys.stderr)

        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            for fut in [ex.submit(one, j) for j in jobs]:
                fut.result()


if __name__ == "__main__":
    main()
