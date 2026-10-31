#!/usr/bin/env python3
"""Wall-clock and peak RSS of C MAFFT, mafft-rs and pymafft from n = 100 to 5000 (issue #22).

For every (set, n, mode, threads) cell of the scaling corpus
(`benchmarks/scaling_corpus.py`) this runs, interleaved and with the
order rotated per replicate:

* `c`     the in-tree C build (`mafft-upstream/scripts/mafft` with
          `MAFFT_BINARIES=mafft-upstream/binaries`), as a subprocess;
* `rust`  the `mafft-rs` CLI (`target/release/mafft-rs`), as a subprocess;
* `py`    native `pymafft.align()` on an in-memory list of
          (name, sequence) tuples, in a fresh interpreter of the venv
          (`.venv/bin/python`). `wall_s` is the time around the `align()`
          call alone; `proc_wall_s` includes interpreter start-up and
          reading the FASTA.

and records wall-clock, peak RSS of the process tree (`ru_maxrss` from
`wait4`, the number `/usr/bin/time -l` / `-v` report), output width,
SHA-256 of the output, and whether the C and Rust outputs are
byte-identical (`identical_c_rust`) and pymafft's rows equal C's
(`identical_c_py`, name and gapped sequence per row; pymafft does not
emit the CLI's line wrapping). Every row also carries the host, `uname
-m`, CPU, FP-contraction policy and the 1-minute load average when the
run started, so TSVs from different platforms can be joined (each
platform's Rust is compared with that platform's C build; `compare`
reports where the two C builds themselves disagree).

Modes (all runs add `--quiet` and `--thread T`):

    fftns2    (no flags)                       FFT-NS-2
    fftnsi    --maxiterate 1000                FFT-NS-i
    linsi     --localpair --maxiterate 1000    L-INS-i
    ginsi     --globalpair --maxiterate 1000   G-INS-i   (n <= 250 only)
    parttree  --parttree                       PartTree  (no pymafft: not exposed by align())

Grid: every set at each of its corpus sizes; sizes above 1000 (5000 for
PF00069 / PF00005, 2000 for 16S) only for fftns2 and parttree.
`--thread 1` everywhere, plus `--thread 8` for n >= 500.

Skip rule. Cells are run in increasing n. Before a cell, the C time is
projected from the largest smaller n of the same (set, mode) already
measured (same thread count if available, else `--thread 1`) as
t * (n / n_prev) ** k with k = 2 (k = 1.2 for parttree). If the projection
exceeds `--budget` (default 3600 s), or C timed out or was skipped at a
smaller n, the cell is recorded as `skipped` for all implementations and
so are the larger n. Each single run has a hard `--timeout` (default
7200 s). In practice this drops L-INS-i / FFT-NS-i at n = 1000 for the
long-sequence families if the n = 500 point already takes > 15 min.

USAGE
-----

    make -C mafft-upstream/core
    cargo build --release -p mafft-rs
    python3.12 -m venv .venv && .venv/bin/pip install maturin biopython
    (cd crates/pymafft && ../../.venv/bin/maturin develop --release)
    benchmarks/scaling_corpus.py build

    # Pipeline check: PF00069, n = 100, FFT-NS-2, --thread 1, 1 replicate
    benchmarks/scaling.py run --smoke

    # Full sweep (resumable: each run is cached under target/scaling/<host>/;
    # re-running fills in what is missing). Run on a quiet machine.
    benchmarks/scaling.py run --reps 3

    # Subsets
    benchmarks/scaling.py run --sets PF00069 16S --sizes 100 250 --modes fftns2 linsi

    # Median + IQR per cell, Rust/C and pymafft/C ratios
    benchmarks/scaling.py summarize benchmarks/results/scaling_<host>_<date>.tsv

    # Agreement of two platforms' C builds (and of their Rust builds)
    benchmarks/scaling.py compare a.tsv b.tsv

`run` writes `benchmarks/results/scaling_<host>_<date>.tsv` (or
`scaling_smoke_<date>.tsv` with `--smoke`), merging with an existing
file of the same name.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import signal
import statistics
import subprocess
import sys
import threading
import time
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "benchmarks"))
import scaling_corpus as corpus  # noqa: E402

MODES = {
    "fftns2": ("FFT-NS-2", [], ("fftns2", 0)),
    "fftnsi": ("FFT-NS-i", ["--maxiterate", "1000"], ("fftnsi", 1000)),
    "linsi": ("L-INS-i", ["--localpair", "--maxiterate", "1000"], ("linsi", 1000)),
    "ginsi": ("G-INS-i", ["--globalpair", "--maxiterate", "1000"], ("ginsi", 1000)),
    "parttree": ("PartTree", ["--parttree"], None),
}
LARGE_N_MODES = ("fftns2", "parttree")
GINSI_MAX_N = 250
MT_MIN_N = 500
GROWTH = {"parttree": 1.2}
IMPLS = ("c", "rust", "py")
COLUMNS = ["host", "uname_m", "cpu", "fp_policy", "set", "molecule", "n", "mode",
           "threads", "impl", "rep", "status", "wall_s", "proc_wall_s", "max_rss_mb",
           "aln_width", "sha256", "identical_c_rust", "identical_c_py", "load1_start"]


# ---------------------------------------------------------------- host info

def cpu_model() -> str:
    try:
        if sys.platform == "darwin":
            return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                  capture_output=True, text=True).stdout.strip()
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "NA"


def default_fp_policy() -> str:
    # mafft-rs mirrors the host C compiler's default: clang/arm64 contracts
    # a*b+c into FMAs, gcc/x86-64 does not.
    return "fma" if platform.machine() in ("arm64", "aarch64") else "none"


# ------------------------------------------------------------------- running

def run_measured(argv: list[str], env: dict, out_path: Path, err_path: Path,
                 timeout: float, result: Path | None = None) -> dict:
    """Run argv with stdout to out_path; wall-clock and ru_maxrss via wait4."""
    with out_path.open("wb") as out, err_path.open("wb") as err:
        t0 = time.perf_counter()
        p = subprocess.Popen(argv, stdout=out, stderr=err, env=env,
                             start_new_session=True)
        timed_out = threading.Event()

        def kill():
            timed_out.set()
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

        timer = threading.Timer(timeout, kill)
        timer.start()
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.perf_counter() - t0
        timer.cancel()
        p.returncode = os.waitstatus_to_exitcode(status)
    # ru_maxrss: bytes on macOS, KiB on Linux.
    rss = ru.ru_maxrss / (1 << 20) if sys.platform == "darwin" else ru.ru_maxrss / 1024
    if timed_out.is_set():
        st = "timeout"
    elif p.returncode != 0:
        tail = " | ".join(err_path.read_text(errors="replace").strip().splitlines()[-2:])
        st = f"exit_{p.returncode}: {tail}"[:200]
    elif not (result or out_path).is_file() or (result or out_path).stat().st_size == 0:
        st = "no_output"
    else:
        st = "ok"
    return {"status": st, "proc_wall_s": round(wall, 4), "max_rss_mb": round(rss, 1)}


def py_child(args) -> int:
    """Child process for the `py` implementation (`scaling.py _py ...`)."""
    import pymafft
    recs = [(h.split()[0], s) for h, s in corpus.parse_fasta(Path(args.fasta).read_text())]
    strategy, maxiterate = MODES[args.mode][2]
    t0 = time.perf_counter()
    res = pymafft.align(recs, strategy=strategy, maxiterate=maxiterate,
                        threads=args.threads, quiet=True)
    el = time.perf_counter() - t0
    Path(args.out).write_text(res.to_fasta())
    print(json.dumps({"wall_s": round(el, 4)}), file=sys.stderr)
    return 0


def rows_of(text: str) -> list[tuple[str, str]]:
    return [(h.split()[0], s.upper()) for h, s in corpus.parse_fasta(text)]


def run_one(args, fasta: Path, set_id: str, n: int, mode: str, threads: int,
            impl: str, rep: int) -> dict:
    d = args.cache / set_id / f"n{n}" / mode / f"t{threads}"
    meta_p = d / f"{impl}.rep{rep}.json"
    if meta_p.is_file():
        meta = json.loads(meta_p.read_text())
        if meta["status"] in ("ok", "timeout", "skipped", "na") or args.no_retry:
            return meta
    d.mkdir(parents=True, exist_ok=True)
    if impl == "py" and MODES[mode][2] is None:
        meta = {"status": "na"}
        meta_p.write_text(json.dumps(meta))
        return meta
    flags = ["--quiet", "--thread", str(threads), *MODES[mode][1]]
    env = dict(os.environ)
    out = d / f"{impl}.rep{rep}.fa"
    err = d / f"{impl}.rep{rep}.err"
    if impl == "c":
        env["MAFFT_BINARIES"] = str(args.mafft_binaries)
        argv = [str(args.mafft), *flags, str(fasta)]
    elif impl == "rust":
        argv = [str(args.mafft_rs), *flags, str(fasta)]
    else:
        argv = [str(args.python), str(Path(__file__).resolve()), "_py", str(fasta),
                str(out), mode, str(threads)]
        out = d / f"{impl}.rep{rep}.stdout"
    load1 = os.getloadavg()[0]
    meta = run_measured(argv, env, out, err, args.timeout,
                        result=d / f"{impl}.rep{rep}.fa")
    meta["load1_start"] = round(load1, 2)
    if impl == "py":
        out.unlink(missing_ok=True)
        out = d / f"{impl}.rep{rep}.fa"
        if meta["status"] == "ok":
            try:
                meta["wall_s"] = json.loads(err.read_text().strip().splitlines()[-1])["wall_s"]
            except (ValueError, IndexError, KeyError):
                meta["status"] = "bad_child_output"
    else:
        meta["wall_s"] = meta["proc_wall_s"]
    if meta["status"] == "ok":
        text = out.read_text()
        meta["sha256"] = hashlib.sha256(text.encode()).hexdigest()
        r = rows_of(text)
        meta["aln_width"] = len(r[0][1]) if r else 0
        meta["rows_sha256"] = hashlib.sha256(
            "\n".join(f"{a}\t{b}" for a, b in r).encode()).hexdigest()
        if rep != 0:
            out.unlink()  # keep one output per cell; hashes cover the others
    meta_p.write_text(json.dumps(meta))
    return meta


def skip_cell(args, set_id: str, n: int, mode: str, threads: int) -> None:
    d = args.cache / set_id / f"n{n}" / mode / f"t{threads}"
    d.mkdir(parents=True, exist_ok=True)
    for impl in IMPLS:
        for rep in range(args.reps):
            p = d / f"{impl}.rep{rep}.json"
            if not p.is_file():
                p.write_text(json.dumps({"status": "skipped"}))


def c_history(args, set_id: str, mode: str, threads: int, below: int,
              sizes: list[int]) -> tuple[int, float | None, bool]:
    """(n_prev, median C wall at n_prev, C failed/skipped at some smaller n)."""
    best: tuple[int, float | None] = (0, None)
    for n in sizes:
        if n >= below:
            continue
        d = args.cache / set_id / f"n{n}" / mode / f"t{threads}"
        metas = [json.loads(p.read_text()) for p in sorted(d.glob("c.rep*.json"))]
        if any(m["status"] in ("timeout", "skipped") for m in metas):
            return n, None, True
        ok = [m["wall_s"] for m in metas if m["status"] == "ok"]
        if ok:
            best = (n, statistics.median(ok))
    return best[0], best[1], False


def should_skip(args, set_id, n, mode, threads, sizes) -> str | None:
    for t in dict.fromkeys((threads, 1)):
        n_prev, t_prev, failed = c_history(args, set_id, mode, t, n, sizes)
        if failed:
            return f"C timed out or skipped at n={n_prev} (--thread {t})"
        if t_prev is not None:
            proj = t_prev * (n / n_prev) ** GROWTH.get(mode, 2.0)
            if proj > args.budget:
                return f"C projected {proj:.0f}s > budget {args.budget}s from n={n_prev}"
            return None
    return None


def grid(args, manifest: list[dict]) -> list[tuple[str, str, int, str, int]]:
    sizes_by_set: dict[str, list[int]] = defaultdict(list)
    mol: dict[str, str] = {}
    for r in manifest:
        sizes_by_set[r["set"]].append(int(r["n"]))
        mol[r["set"]] = r["molecule"]
    cells = []
    for set_id in args.sets:
        for n in sorted(sizes_by_set.get(set_id, [])):
            if args.sizes and n not in args.sizes:
                continue
            for mode in args.modes:
                if n > 1000 and mode not in LARGE_N_MODES:
                    continue
                if mode == "ginsi" and n > GINSI_MAX_N:
                    continue
                for t in args.threads:
                    if t != 1 and n < MT_MIN_N:
                        continue
                    cells.append((set_id, mol[set_id], n, mode, t))
    order = {m: i for i, m in enumerate(MODES)}
    cells.sort(key=lambda c: (c[2], order[c[3]], corpus.ALL_SETS.index(c[0]), c[4]))
    return cells


def cmd_run(args) -> int:
    if args.smoke:
        args.sets, args.sizes, args.modes = ["PF00069"], [100], ["fftns2"]
        args.threads, args.reps = [1], 1
    for p, what in ((args.mafft, "C mafft wrapper (make -C mafft-upstream/core)"),
                    (args.mafft_rs, "mafft-rs (cargo build --release -p mafft-rs)"),
                    (args.python, "venv python with pymafft (see USAGE)")):
        if not Path(p).is_file():
            print(f"error: {what} not found at {p}", file=sys.stderr)
            return 2
    manifest = corpus.read_manifest()
    host = args.host
    args.cache = args.cache / host
    info = {"host": host, "uname_m": platform.machine(), "cpu": cpu_model(),
            "fp_policy": args.fp_policy}
    cells = grid(args, manifest)
    missing = {c[0] for c in cells
               if not (corpus.DATA_DIR / c[0] / f"n{c[2]}.fa").is_file()}
    if missing:
        print(f"error: corpus files missing for {sorted(missing)}; run "
              "benchmarks/scaling_corpus.py build", file=sys.stderr)
        return 2
    print(f"{info}; load {os.getloadavg()}; {len(cells)} cells x {len(IMPLS)} impls x "
          f"{args.reps} reps", file=sys.stderr)
    sizes_by_set = defaultdict(list)
    for r in manifest:
        sizes_by_set[r["set"]].append(int(r["n"]))

    t_start = time.time()
    for k, (set_id, mol, n, mode, threads) in enumerate(cells, 1):
        fasta = corpus.DATA_DIR / set_id / f"n{n}.fa"
        why = should_skip(args, set_id, n, mode, threads, sizes_by_set[set_id])
        if why:
            skip_cell(args, set_id, n, mode, threads)
            print(f"[{k}/{len(cells)}] {set_id} n={n} {mode} t{threads}: skipped ({why})",
                  file=sys.stderr)
            continue
        for rep in range(args.reps):
            order = IMPLS[rep % 3:] + IMPLS[:rep % 3]
            for impl in order:
                m = run_one(args, fasta, set_id, n, mode, threads, impl, rep)
                if not args.quiet:
                    print(f"[{k}/{len(cells)} {time.time() - t_start:7.0f}s] {set_id} n={n} "
                          f"{mode} t{threads} rep{rep} {impl:<4} {m['status']} "
                          f"{m.get('wall_s', 0):.2f}s {m.get('max_rss_mb', 0):.0f}MB",
                          file=sys.stderr)

    rows = collect(args, cells, info)
    out = args.output or (REPO / "benchmarks" / "results" /
                          (f"scaling_smoke_{dt.date.today().isoformat()}.tsv" if args.smoke
                           else f"scaling_{host}_{dt.date.today().isoformat()}.tsv"))
    write_tsv(out, rows)
    print(f"wrote {out}", file=sys.stderr)
    bad = [r for r in rows if r["impl"] == "rust" and r["identical_c_rust"] == "no"]
    print(summarize(rows))
    if bad:
        print(f"{len(bad)} Rust outputs differ from C", file=sys.stderr)
    return 0


def collect(args, cells, info) -> list[dict]:
    rows = []
    for set_id, mol, n, mode, threads in cells:
        d = args.cache / set_id / f"n{n}" / mode / f"t{threads}"
        for rep in range(args.reps):
            metas = {}
            for impl in IMPLS:
                p = d / f"{impl}.rep{rep}.json"
                metas[impl] = json.loads(p.read_text()) if p.is_file() else {"status": "missing"}
            c = metas["c"]
            for impl in IMPLS:
                m = metas[impl]
                row = {**info, "set": set_id, "molecule": mol, "n": n, "mode": mode,
                       "threads": threads, "impl": impl, "rep": rep, "status": m["status"],
                       "wall_s": m.get("wall_s", "NA"), "proc_wall_s": m.get("proc_wall_s", "NA"),
                       "max_rss_mb": m.get("max_rss_mb", "NA"),
                       "aln_width": m.get("aln_width", "NA"), "sha256": m.get("sha256", "NA"),
                       "identical_c_rust": "NA", "identical_c_py": "NA",
                       "load1_start": m.get("load1_start", "NA")}
                r_m, p_m = metas["rust"], metas["py"]
                if c["status"] == "ok" and r_m["status"] == "ok":
                    row["identical_c_rust"] = "yes" if c["sha256"] == r_m["sha256"] else "no"
                if c["status"] == "ok" and p_m["status"] == "ok":
                    row["identical_c_py"] = ("yes" if c["rows_sha256"] == p_m["rows_sha256"]
                                             else "no")
                rows.append(row)
    return rows


def key(r: dict) -> tuple:
    return (r["host"], r["set"], str(r["n"]), r["mode"], str(r["threads"]), r["impl"],
            str(r["rep"]))


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    merged = {key(r): r for r in (read_tsv(path) if path.is_file() else [])}
    merged.update({key(r): {k: str(v) for k, v in r.items()} for r in rows})
    order = {m: i for i, m in enumerate(MODES)}
    sets = {s: i for i, s in enumerate(corpus.ALL_SETS)}
    with path.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for k in sorted(merged, key=lambda k: (k[0], int(k[2]), order.get(k[3], 9),
                                               sets.get(k[1], 99), int(k[4]),
                                               IMPLS.index(k[5]), int(k[6]))):
            fh.write("\t".join(merged[k][c] for c in COLUMNS) + "\n")


# ----------------------------------------------------------------- reporting

def quartiles(xs: list[float]) -> tuple[float, float, float]:
    xs = sorted(xs)
    if len(xs) == 1:
        return xs[0], xs[0], xs[0]
    q = statistics.quantiles(xs, n=4, method="inclusive")
    return q[1], q[0], q[2]


def summarize(rows: list[dict]) -> str:
    cells: dict[tuple, dict] = defaultdict(lambda: defaultdict(list))
    ident: dict[tuple, set] = defaultdict(set)
    rss: dict[tuple, dict] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        k = (r["set"], int(r["n"]), r["mode"], int(r["threads"]))
        if r["status"] == "ok":
            cells[k][r["impl"]].append(float(r["wall_s"]))
            rss[k][r["impl"]].append(float(r["max_rss_mb"]))
        elif r["impl"] == "c":
            cells[k]["_status"].append(r["status"])
        if r["impl"] == "rust":
            ident[k].add(r["identical_c_rust"])
    lines = ["| Set | n | Mode | T | C s (median [IQR]) | Rust s | pymafft s | Rust/C "
             "| py/C | RSS C / Rust MB | C = Rust |",
             "|---|---:|---|---:|---|---|---|---:|---:|---|---|"]
    order = {m: i for i, m in enumerate(MODES)}
    for k in sorted(cells, key=lambda k: (corpus.ALL_SETS.index(k[0])
                                          if k[0] in corpus.ALL_SETS else 99,
                                          k[1], order[k[2]], k[3])):
        c = cells[k]

        def fmt(impl):
            if not c.get(impl):
                return "–"
            m, lo, hi = quartiles(c[impl])
            return f"{m:.2f} [{lo:.2f}–{hi:.2f}]"

        def ratio(impl):
            if not c.get(impl) or not c.get("c"):
                return "–"
            return f"{statistics.median(c[impl]) / statistics.median(c['c']):.2f}"

        if not c.get("c") and c.get("_status"):
            lines.append(f"| {k[0]} | {k[1]} | {MODES[k[2]][0]} | {k[3]} | "
                         f"{c['_status'][0]} | | | | | | |")
            continue
        r_c = statistics.median(rss[k]["c"]) if rss[k]["c"] else float("nan")
        r_r = statistics.median(rss[k]["rust"]) if rss[k]["rust"] else float("nan")
        idv = ident[k] - {"NA"}
        lines.append(f"| {k[0]} | {k[1]} | {MODES[k[2]][0]} | {k[3]} | {fmt('c')} | "
                     f"{fmt('rust')} | {fmt('py')} | {ratio('rust')} | {ratio('py')} | "
                     f"{r_c:.0f} / {r_r:.0f} | "
                     f"{'yes' if idv == {'yes'} else ('no' if 'no' in idv else '–')} |")
    return "\n".join(lines)


def cmd_summarize(args) -> int:
    print(summarize(read_tsv(Path(args.tsv))))
    return 0


def cmd_compare(args) -> int:
    """Join two hosts' TSVs on (set, n, mode, threads, rep 0)."""
    def index(path):
        out = {}
        for r in read_tsv(Path(path)):
            if r["rep"] == "0" and r["status"] == "ok":
                out[(r["set"], r["n"], r["mode"], r["threads"], r["impl"])] = r
        return out

    a, b = index(args.a), index(args.b)
    ha = next(iter(a.values()))["host"] if a else "?"
    hb = next(iter(b.values()))["host"] if b else "?"
    print(f"A = {args.a} ({ha}), B = {args.b} ({hb})")
    print("| Set | n | Mode | T | C(A) = C(B) | width A / B | Rust(A) = C(A) | Rust(B) = C(B) |")
    print("|---|---:|---|---:|---|---|---|---|")
    same = total = 0
    for k in sorted(set(a) & set(b)):
        if k[4] != "c":
            continue
        ra, rb = a.get(k[:4] + ("rust",)), b.get(k[:4] + ("rust",))
        eq = a[k]["sha256"] == b[k]["sha256"]
        total += 1
        same += eq
        print(f"| {k[0]} | {k[1]} | {k[2]} | {k[3]} | {'yes' if eq else 'no'} | "
              f"{a[k]['aln_width']} / {b[k]['aln_width']} | "
              f"{ra['identical_c_rust'] if ra else '–'} | {rb['identical_c_rust'] if rb else '–'} |")
    print(f"\nC builds agree on {same} of {total} common cells")
    return 0


# ------------------------------------------------------------------------ CLI

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("Modes")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="time the grid (cached, resumable)")
    r.add_argument("--smoke", action="store_true",
                   help="PF00069, n=100, FFT-NS-2, --thread 1, one replicate")
    r.add_argument("--sets", nargs="+", default=corpus.ALL_SETS, choices=corpus.ALL_SETS)
    r.add_argument("--sizes", nargs="+", type=int, default=None)
    r.add_argument("--modes", nargs="+", default=list(MODES), choices=list(MODES))
    r.add_argument("--threads", nargs="+", type=int, default=[1, 8])
    r.add_argument("--reps", type=int, default=3)
    r.add_argument("--timeout", type=float, default=7200, help="seconds per run")
    r.add_argument("--budget", type=float, default=3600,
                   help="skip a cell if C's projected time exceeds this (seconds)")
    r.add_argument("--host", default=f"{platform.system().lower()}-{platform.machine()}",
                   help="label for the cache directory and the results file")
    r.add_argument("--fp-policy", default=default_fp_policy(),
                   help="FP-contraction policy of both builds (default: host default)")
    r.add_argument("--mafft", type=Path, default=REPO / "mafft-upstream/scripts/mafft")
    r.add_argument("--mafft-binaries", type=Path, default=REPO / "mafft-upstream/binaries")
    r.add_argument("--mafft-rs", type=Path, default=REPO / "target/release/mafft-rs")
    r.add_argument("--python", type=Path, default=REPO / ".venv/bin/python",
                   help="interpreter with pymafft installed")
    r.add_argument("--cache", type=Path, default=REPO / "target/scaling")
    r.add_argument("--output", type=Path, default=None)
    r.add_argument("--no-retry", action="store_true",
                   help="do not retry cached failures")
    r.add_argument("--quiet", action="store_true")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("summarize", help="Markdown table of a results TSV")
    s.add_argument("tsv")
    s.set_defaults(func=cmd_summarize)

    c = sub.add_parser("compare", help="C-vs-C agreement of two platforms' TSVs")
    c.add_argument("a")
    c.add_argument("b")
    c.set_defaults(func=cmd_compare)

    y = sub.add_parser("_py", help=argparse.SUPPRESS)
    y.add_argument("fasta")
    y.add_argument("out")
    y.add_argument("mode", choices=list(MODES))
    y.add_argument("threads", type=int)
    y.set_defaults(func=py_child)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
