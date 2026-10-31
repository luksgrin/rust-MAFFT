#!/usr/bin/env python3
"""BAliBASE 3 SP / TC scores of C MAFFT and mafft-rs (issue #28).

Byte-identity is the correctness criterion for rust-MAFFT; sum-of-pairs
(SP) and total-column (TC) scores against curated references are how
alignment *methods* are usually compared. This script computes the
latter for both implementations so the two views can be put side by
side: for every BAliBASE 3 set and mode it aligns the `.tfa` input with
the in-tree C MAFFT and with `mafft-rs`, scores each output against the
set's `.xml` reference (core blocks only), and records whether the two
outputs are byte-identical.

Scoring uses BAliBASE's own `bali_score` (v3.01), compiled from
`$BALIBASE_DIR/bali_score_src` with the system `cc` and expat on first
use. A Python re-implementation of the same definitions (`py_score`) runs
alongside it on every alignment and the run fails loudly if the two
disagree, so the scorer is validated on the whole corpus rather than on
a handful of sets.

Definitions (bali_score, XML reference, core blocks):

* SP = (residue pairs aligned in the test alignment that are also
  aligned in a core-block column of the reference) / (number of pairs in
  the reference core-block columns, counted as N*(N-1)/2 per column).
* TC = (test columns in which all N sequences carry the residues of one
  and the same reference core-block column) / (test columns in which the
  first reference sequence carries a core-block residue), truncated to
  two decimals as bali_score does.

USAGE
-----

    make -C mafft-upstream/core
    cargo build --release -p mafft-rs
    export BALIBASE_DIR=/path/to/bb3_release      # contains RV11/ ... RV50/

    # Align and score (resumable: alignments are cached under
    # target/balibase_scores/, re-running only fills in what is missing):
    benchmarks/balibase_scores.py run --modes fftns2
    benchmarks/balibase_scores.py run --modes fftnsi linsi ginsi einsi --jobs 6
    benchmarks/balibase_scores.py run --filter RV11 --modes linsi

    # Per-(subset, mode) summary of a results file, as a Markdown table:
    benchmarks/balibase_scores.py summarize benchmarks/results/balibase_scores_2026-10-31.tsv

    # The two counter-examples on mafft-upstream/test/sample (needs a
    # second mafft-rs built with --features fp-contract-none):
    CARGO_TARGET_DIR=target/nofma cargo build --release -p mafft-rs --features fp-contract-none
    benchmarks/balibase_scores.py counterexamples

`run` writes `benchmarks/results/balibase_scores_<date>.tsv` with columns
`set subset mode impl n ref_len aln_len SP TC bytes_identical`
(`bytes_identical` compares the C and Rust outputs of the same set and
mode; SP/TC are `NA` if an implementation failed or timed out).

Modes (the flag sets are the ones used by the byte-identity sweep):

    fftns2  (no flags)                       FFT-NS-2
    fftnsi  --maxiterate 1000                FFT-NS-i
    linsi   --localpair --maxiterate 1000    L-INS-i
    ginsi   --globalpair --maxiterate 1000   G-INS-i
    einsi   --genafpair --maxiterate 1000    E-INS-i
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

MODES = {
    "fftns2": ("FFT-NS-2", []),
    "fftnsi": ("FFT-NS-i", ["--maxiterate", "1000"]),
    "linsi": ("L-INS-i", ["--localpair", "--maxiterate", "1000"]),
    "ginsi": ("G-INS-i", ["--globalpair", "--maxiterate", "1000"]),
    "einsi": ("E-INS-i", ["--genafpair", "--maxiterate", "1000"]),
}
IMPLS = ("c", "rust")
SUBSETS = ("RV11", "RV12", "RV20", "RV30", "RV40", "RV50")
COLUMNS = ["set", "subset", "mode", "impl", "n", "ref_len", "aln_len",
           "SP", "TC", "bytes_identical"]


# ---------------------------------------------------------------- FASTA / MSF

def read_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA; names are the first whitespace-separated header token."""
    out: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            out.append((line[1:].split()[0] if line[1:].split() else "", []))
        elif out:
            out[-1][1].append(line.strip())
    return [(n, "".join(s)) for n, s in out]


def write_msf(seqs: list[tuple[str, str]], path: Path) -> None:
    """Minimal GCG MSF that bali_score's `read_msf` accepts."""
    width = max(len(s) for _, s in seqs)
    lines = [f" MSF: {width}  Type: P  Check: 0 ..", ""]
    for name, s in seqs:
        lines.append(f" Name: {name} Len: {width}  Check: 0  Weight: 1.0")
    lines += ["", "//", ""]
    for start in range(0, width, 50):
        for name, s in seqs:
            lines.append(f"{name:<32}{s[start:start + 50].replace('-', '.')}")
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------ reference (XML)

class Reference:
    """BAliBASE XML reference: aligned sequences plus the core-block mask."""

    def __init__(self, xml_path: Path):
        root = ET.parse(xml_path).getroot()
        aln = root.find("alignment")
        self.names: list[str] = []
        self.rows: list[str] = []
        for seq in aln.findall("sequence"):
            self.names.append(seq.findtext("seq-name").strip())
            data = re.sub(r"\s+", "", seq.findtext("seq-data"))
            self.rows.append(data.replace(".", "-").replace("~", "-"))
        self.core: list[bool] | None = None
        for cs in aln.iter("column-score"):
            if cs.findtext("colsco-name", "").strip() == "coreblock":
                vals = [int(v) for v in cs.findtext("colsco-data").split()]
                self.core = [v == 1 for v in vals]
                break
        self.length = max(len(r) for r in self.rows)


def py_score(ref: Reference, test: list[tuple[str, str]]) -> tuple[float, float]:
    """SP and TC with bali_score 3.01's definitions (core-block mode)."""
    n = len(ref.names)
    by_name = {name.lower(): s for name, s in test}
    if ref.core is None:
        raise ValueError("reference has no coreblock annotation")
    core = ref.core + [False] * max(0, ref.length - len(ref.core))
    # Code each residue of each test row (reference order) with its
    # 1-based reference column if that column is a core block, else 0.
    coded: list[list[int]] = []
    for name, refrow in zip(ref.names, ref.rows):
        cols = [i for i, c in enumerate(refrow) if c != "-"]
        row = by_name[name.lower()]
        out, k = [], 0
        for c in row:
            if c == "-":
                out.append(0)
            else:
                col = cols[k]
                out.append(col + 1 if core[col] else 0)
                k += 1
        coded.append(out)
    width = max(len(r) for r in coded)
    for r in coded:
        r.extend([0] * (width - len(r)))
    max_pairs = sum(1 for c in core[:ref.length] if c) * (n * (n - 1) // 2)
    sp = tc = ncols = 0
    for i in range(width):
        counts: dict[int, int] = defaultdict(int)
        for r in coded:
            if r[i]:
                counts[r[i]] += 1
        sp += sum(c * (c - 1) // 2 for c in counts.values())
        first = coded[0][i]
        if first:
            ncols += 1
            if counts[first] >= n:
                tc += 1
    tc_pct = int(100 * tc / ncols) if ncols else 0
    return sp / max_pairs, tc_pct / 100.0


# ----------------------------------------------------------------- bali_score

def ensure_bali_score(balibase: Path, cache: Path) -> Path | None:
    exe = cache / "bali_score"
    if exe.is_file():
        return exe
    src = balibase / "bali_score_src"
    if not (src / "bali_score.c").is_file():
        return None
    build = cache / "bali_score_build"
    build.mkdir(parents=True, exist_ok=True)
    for f in src.glob("*.[ch]"):
        shutil.copy(f, build / f.name)
    cfiles = ["readxml.c", "init.c", "util.c", "bali_score.c"]
    # macOS enables _FORTIFY_SOURCE by default and its object-size check
    # traps on bali_score's strcpy into struct fields (the code is clean
    # under -fsanitize=address), so turn fortification off.
    subprocess.run(["cc", "-O", "-w", "-U_FORTIFY_SOURCE", "-D_FORTIFY_SOURCE=0",
                    "-o", str(exe), *cfiles, "-lexpat", "-lm"],
                   cwd=build, check=True)
    return exe


def bali_score(exe: Path, xml: Path, test: list[tuple[str, str]]) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        msf = Path(td) / "test.msf"
        write_msf(test, msf)
        proc = subprocess.run([str(exe), str(xml), str(msf)],
                              capture_output=True, text=True, check=False)
    m = re.search(r"^auto \S+ (\S+) (\S+)$", proc.stdout, re.M)
    if not m:
        raise RuntimeError(f"bali_score failed on {xml}: {proc.stdout[-500:]}"
                           f"{proc.stderr[-500:]}")
    return float(m.group(1)), float(m.group(2))


# ------------------------------------------------------------------- aligning

def run_aligner(argv: list[str], env: dict, timeout: int) -> tuple[str | None, float, str]:
    t0 = time.perf_counter()
    try:
        p = subprocess.run(argv, capture_output=True, text=True, env=env,
                           timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return None, time.perf_counter() - t0, "timeout"
    el = time.perf_counter() - t0
    if p.returncode != 0:
        tail = " | ".join(p.stderr.strip().splitlines()[-3:])
        return None, el, f"exit_{p.returncode}: {tail}"
    if not p.stdout:
        return None, el, "no_output"
    return p.stdout, el, ""


def align_cached(args, tfa: Path, mode: str, impl: str) -> tuple[str | None, dict]:
    """Return (aligned FASTA or None, metadata); reuse the cache if present."""
    d = args.cache / mode / impl
    out = d / f"{tfa.stem}.fa"
    meta_p = d / f"{tfa.stem}.json"
    if meta_p.is_file():
        meta = json.loads(meta_p.read_text())
        # Successes and timeouts are final; other failures (a non-zero
        # exit, e.g. from a transient resource problem) are retried.
        if meta["status"] == "ok":
            return out.read_text(), meta
        if meta["status"] == "timeout" or args.no_retry:
            return None, meta
    flags = MODES[mode][1]
    env = dict(os.environ)
    if impl == "c":
        env["MAFFT_BINARIES"] = str(args.mafft_binaries)
        argv = [str(args.mafft), *flags, str(tfa)]
    elif impl == "rust-nofma":
        argv = [str(args.mafft_rs_nofma), *flags, str(tfa)]
    else:
        argv = [str(args.mafft_rs), *flags, str(tfa)]
    text, el, err = run_aligner(argv, env, args.timeout)
    meta = {"status": "ok" if text is not None else err, "elapsed": round(el, 3),
            "argv": " ".join(argv[1:-1])}
    d.mkdir(parents=True, exist_ok=True)
    if text is not None:
        out.write_text(text)
    meta_p.write_text(json.dumps(meta))
    return text, meta


def discover(balibase: Path, filt: str) -> list[Path]:
    tests = []
    for sub in SUBSETS:
        tests += sorted((balibase / sub).glob("BB*.tfa"))
    return [t for t in tests if filt in str(t)]


def cmd_run(args) -> int:
    if not args.balibase or not Path(args.balibase).is_dir():
        print("error: set BALIBASE_DIR or pass --balibase", file=sys.stderr)
        return 2
    balibase = Path(args.balibase)
    for p, what in ((args.mafft, "C mafft wrapper (make -C mafft-upstream/core)"),
                    (args.mafft_rs, "mafft-rs (cargo build --release -p mafft-rs)")):
        if not Path(p).is_file():
            print(f"error: {what} not found at {p}", file=sys.stderr)
            return 2
    args.cache.mkdir(parents=True, exist_ok=True)
    exe = ensure_bali_score(balibase, args.cache)
    print(f"scorer: {'bali_score ' + str(exe) if exe else 'Python only (no bali_score_src)'}",
          file=sys.stderr)
    tests = discover(balibase, args.filter)
    jobs = [(t, m, i) for m in args.modes for t in tests for i in IMPLS]
    print(f"{len(tests)} sets x {len(args.modes)} modes x 2 impls = {len(jobs)} alignments",
          file=sys.stderr)

    done = 0
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(align_cached, args, t, m, i): (t, m, i) for t, m, i in jobs}
        for f in cf.as_completed(futs):
            t, m, i = futs[f]
            _, meta = f.result()
            done += 1
            if not args.quiet:
                print(f"[{done:5d}/{len(jobs)} {time.time() - t0:7.0f}s] {m:<7} {i:<4} "
                      f"{t.stem} {meta['status']} {meta['elapsed']:.2f}s", file=sys.stderr)

    # Score everything that is cached for the requested scope (no second
    # attempt at failures here; that is the next invocation's job).
    args.no_retry = True
    rows = []
    scorer_mismatch = []
    refs: dict[Path, Reference] = {}
    for m in args.modes:
        for t in tests:
            xml = t.with_suffix(".xml")
            ref = refs.setdefault(xml, Reference(xml))
            outs = {}
            for i in IMPLS:
                outs[i], _ = align_cached(args, t, m, i)
            ident = (outs["c"] is not None and outs["c"] == outs["rust"])
            for i in IMPLS:
                row = {"set": t.stem, "subset": t.parent.name, "mode": m, "impl": i,
                       "n": len(ref.names), "ref_len": ref.length}
                if outs[i] is None:
                    row.update(aln_len="NA", SP="NA", TC="NA")
                else:
                    aln = read_fasta(outs[i])
                    psp, ptc = py_score(ref, aln)
                    if exe:
                        sp, tc = bali_score(exe, xml, aln)
                        if abs(sp - psp) > 0.0011 or abs(tc - ptc) > 1e-9:
                            scorer_mismatch.append((t.stem, m, i, sp, psp, tc, ptc))
                    else:
                        sp, tc = round(psp, 3), ptc
                    row.update(aln_len=len(aln[0][1]), SP=f"{sp:.3f}", TC=f"{tc:.3f}")
                row["bytes_identical"] = "yes" if ident else "no"
                rows.append(row)

    out = args.output or (REPO / "benchmarks" / "results" /
                          f"balibase_scores_{dt.date.today().isoformat()}.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    # Merge with an existing file so separate invocations (per mode or
    # per subset) accumulate into one table.
    merged: dict[tuple, dict] = {}
    if out.is_file() and not args.overwrite:
        with out.open() as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                merged[(r["mode"], r["set"], r["impl"])] = r
    for r in rows:
        merged[(r["mode"], r["set"], r["impl"])] = {k: str(v) for k, v in r.items()}
    order = {m: k for k, m in enumerate(MODES)}
    with out.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for key in sorted(merged, key=lambda k: (order[k[0]], k[1], k[2])):
            fh.write("\t".join(merged[key][c] for c in COLUMNS) + "\n")
    print(f"wrote {out} ({len(merged)} rows)", file=sys.stderr)

    if scorer_mismatch:
        print("bali_score and py_score DISAGREE on:", file=sys.stderr)
        for x in scorer_mismatch:
            print("  ", *x, file=sys.stderr)
        return 1
    if exe:
        print(f"bali_score and py_score agree on all {sum(r['SP'] != 'NA' for r in rows)} "
              "scored alignments", file=sys.stderr)
    print(summarize(list(merged.values())))
    return 0


# ------------------------------------------------------------------ summarize

def summarize(rows: list[dict]) -> str:
    groups: dict[tuple, dict] = defaultdict(lambda: {"c": [], "rust": [], "ident": 0,
                                                      "differ": 0, "sets": set(), "na": 0})
    for r in rows:
        g = groups[(r["subset"], r["mode"])]
        g["sets"].add(r["set"])
        if r["SP"] == "NA":
            g["na"] += 1
            continue
        g[r["impl"]].append((r["set"], float(r["SP"]), float(r["TC"])))
    for g in groups.values():
        c = {s: (sp, tc) for s, sp, tc in g["c"]}
        rs = {s: (sp, tc) for s, sp, tc in g["rust"]}
        g["differ"] = sum(1 for s in c if s in rs and c[s] != rs[s])
    ident: dict[tuple, set] = defaultdict(set)
    for r in rows:
        if r["bytes_identical"] == "yes":
            ident[(r["subset"], r["mode"])].add(r["set"])

    def mean(xs, k):
        return sum(x[k] for x in xs) / len(xs) if xs else float("nan")

    lines = ["| Subset | Mode | Sets | Byte-identical | SP C | SP Rust | TC C | TC Rust "
             "| SP/TC differ |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    order = {m: k for k, m in enumerate(MODES)}
    tot: dict[str, dict] = defaultdict(lambda: {"c": [], "rust": [], "sets": 0,
                                                "ident": 0, "differ": 0})
    for (sub, m) in sorted(groups, key=lambda k: (order[k[1]], k[0])):
        g = groups[(sub, m)]
        nid = len(ident[(sub, m)])
        lines.append(f"| {sub} | {MODES[m][0]} | {len(g['sets'])} | {nid} | "
                     f"{mean(g['c'], 1):.3f} | {mean(g['rust'], 1):.3f} | "
                     f"{mean(g['c'], 2):.3f} | {mean(g['rust'], 2):.3f} | {g['differ']} |")
        t = tot[m]
        t["c"] += g["c"]
        t["rust"] += g["rust"]
        t["sets"] += len(g["sets"])
        t["ident"] += nid
        t["differ"] += g["differ"]
    for m in sorted(tot, key=order.get):
        t = tot[m]
        lines.append(f"| **all** | **{MODES[m][0]}** | **{t['sets']}** | **{t['ident']}** | "
                     f"**{mean(t['c'], 1):.3f}** | **{mean(t['rust'], 1):.3f}** | "
                     f"**{mean(t['c'], 2):.3f}** | **{mean(t['rust'], 2):.3f}** | "
                     f"**{t['differ']}** |")
    return "\n".join(lines)


def cmd_summarize(args) -> int:
    with open(args.tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh]
    print(summarize(rows))
    return 0


# ------------------------------------------------------------------ crossarch

def cmd_crossarch(args) -> int:
    """Host C build vs the other FP-contraction policy, scored on BAliBASE.

    On arm64 the in-tree C contracts a*b+c into FMAs; `mafft-rs` built with
    `--features fp-contract-none` reproduces the x86-64/gcc C build. Both
    are legitimate builds of the same C source, so wherever their outputs
    differ this measures how far apart SP/TC put two alignments that byte
    comparison rejects.
    """
    balibase = Path(args.balibase)
    if not Path(args.mafft_rs_nofma).is_file():
        print(f"error: {args.mafft_rs_nofma} not found (CARGO_TARGET_DIR=target/nofma "
              "cargo build --release -p mafft-rs --features fp-contract-none)", file=sys.stderr)
        return 2
    exe = ensure_bali_score(balibase, args.cache)
    tests = discover(balibase, args.filter)
    jobs = [(t, m, i) for m in args.modes for t in tests for i in ("c", "rust-nofma")]
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        list(ex.map(lambda j: align_cached(args, *j), jobs))
    args.no_retry = True
    cols = ["set", "subset", "mode", "n", "aln_len_host", "aln_len_other", "SP_host",
            "SP_other", "TC_host", "TC_other", "bytes_identical"]
    rows = []
    for m in args.modes:
        for t in tests:
            xml = t.with_suffix(".xml")
            ref = Reference(xml)
            a, _ = align_cached(args, t, m, "c")
            b, _ = align_cached(args, t, m, "rust-nofma")
            if a is None or b is None:
                continue
            sa, sb = read_fasta(a), read_fasta(b)
            spa, tca = bali_score(exe, xml, sa) if exe else py_score(ref, sa)
            spb, tcb = bali_score(exe, xml, sb) if exe else py_score(ref, sb)
            rows.append([t.stem, t.parent.name, m, len(ref.names), len(sa[0][1]),
                         len(sb[0][1]), f"{spa:.3f}", f"{spb:.3f}", f"{tca:.3f}",
                         f"{tcb:.3f}", "yes" if a == b else "no"])
    out = args.output or (REPO / "benchmarks" / "results" /
                          f"balibase_scores_crossarch_{dt.date.today().isoformat()}.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {out}", file=sys.stderr)
    for m in args.modes:
        rs = [r for r in rows if r[2] == m]
        diff = [r for r in rs if r[10] == "no"]
        print(f"\n{MODES[m][0]}: {len(rs)} sets, {len(diff)} not byte-identical "
              f"between the two C policies")
        if diff:
            print("| Set | Width (host / other) | SP (host / other) | TC (host / other) |")
            print("|---|---|---|---|")
            for r in diff:
                print(f"| {r[0]} | {r[4]} / {r[5]} | {r[6]} / {r[7]} | {r[8]} / {r[9]} |")
            dsp = [abs(float(r[6]) - float(r[7])) for r in diff]
            dtc = [abs(float(r[8]) - float(r[9])) for r in diff]
            same = sum(1 for r in diff if r[6] == r[7] and r[8] == r[9])
            print(f"\nmean |dSP| = {sum(dsp) / len(dsp):.4f} (max {max(dsp):.3f}), "
                  f"mean |dTC| = {sum(dtc) / len(dtc):.4f} (max {max(dtc):.3f}); "
                  f"{same} of {len(diff)} have identical SP and TC")
    return 0


# ------------------------------------------------------------ counterexamples

def residue_columns(aln: list[tuple[str, str]]) -> list[tuple]:
    """Each column as a tuple of residue indices (None for a gap)."""
    counters = [0] * len(aln)
    cols = []
    for i in range(len(aln[0][1])):
        col = []
        for k, (_, s) in enumerate(aln):
            if s[i] == "-":
                col.append(None)
            else:
                col.append(counters[k])
                counters[k] += 1
        cols.append(tuple(col))
    return cols


def aligned_pairs(cols: list[tuple]) -> set:
    pairs = set()
    for col in cols:
        res = [(k, r) for k, r in enumerate(col) if r is not None]
        for a in range(len(res)):
            for b in range(a + 1, len(res)):
                pairs.add((res[a], res[b]))
    return pairs


def compare(a_txt: str, b_txt: str, score_bin: Path) -> dict:
    a, b = read_fasta(a_txt), read_fasta(b_txt)
    ca, cb = residue_columns(a), residue_columns(b)
    pa, pb = aligned_pairs(ca), aligned_pairs(cb)
    shared_cols = len(set(ca) & set(cb))

    def mafft_score(txt):
        # MAFFT's own `score` utility (BLOSUM62, default gap penalties);
        # it rewinds stdin, so it needs a real file.
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as fh:
            fh.write(txt)
            fh.flush()
            with open(fh.name) as inp:
                p = subprocess.run([str(score_bin)], stdin=inp, capture_output=True,
                                   text=True, check=False)
        m = re.search(r"score = (\S+)", p.stdout)
        return float(m.group(1)) if m else float("nan")

    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "a").write_text(a_txt)
        (Path(td) / "b").write_text(b_txt)
        d = subprocess.run(["diff", str(Path(td) / "a"), str(Path(td) / "b")],
                           capture_output=True, text=True, check=False).stdout
    changed = sum(1 for line in d.splitlines() if line[:1] in "<>")
    return {
        "width": (len(ca), len(cb)),
        "identical_columns": shared_cols,
        "tc_a_vs_b": shared_cols / len(cb),
        "tc_b_vs_a": shared_cols / len(ca),
        "pairs": (len(pa), len(pb)),
        "shared_pairs": len(pa & pb),
        "sp_a_vs_b": len(pa & pb) / len(pb),
        "sp_b_vs_a": len(pa & pb) / len(pa),
        "mafft_score_blosum62": (mafft_score(a_txt), mafft_score(b_txt)),
        "changed_lines": changed,
        "bytes_identical": a_txt == b_txt,
    }


def cmd_counterexamples(args) -> int:
    sample = REPO / "mafft-upstream" / "test" / "sample"
    env = dict(os.environ, MAFFT_BINARIES=str(args.mafft_binaries))
    score_bin = args.mafft_binaries / "score"
    cases = [
        ("FMA vs no-FMA (--bl 50 --retree 2 --maxiterate 0)",
         ["--quiet", "--bl", "50", "--retree", "2", "--maxiterate", "0"],
         ("C, arm64 in-tree (fused)", [str(args.mafft)]),
         ("mafft-rs fp-contract-none (= x86-64 C)", [str(args.mafft_rs_nofma)])),
        ("tied-trace residual (--exp 4.30, FFT-NS-2)",
         ["--quiet", "--exp", "4.30"],
         ("C, arm64 in-tree", [str(args.mafft)]),
         ("mafft-rs (default policy)", [str(args.mafft_rs)])),
    ]
    def scoreout(cmd, flags):
        # The run's own objective, as reported by --scoreout (same matrix
        # and penalties as the alignment itself).
        flags = [f for f in flags if f != "--quiet"]
        p = subprocess.run([*cmd, "--scoreout", *flags, str(sample)], capture_output=True,
                           text=True, env=env, check=True)
        m = re.search(r"sum-of-pairs score = (\S+)", p.stderr)
        return float(m.group(1)) if m else float("nan")

    for title, flags, (la, ca), (lb, cbin) in cases:
        a = subprocess.run([*ca, *flags, str(sample)], capture_output=True, text=True,
                           env=env, check=True).stdout
        b = subprocess.run([*cbin, *flags, str(sample)], capture_output=True, text=True,
                           env=env, check=True).stdout
        r = compare(a, b, score_bin)
        r["scoreout"] = (scoreout(ca, flags), scoreout(cbin, flags))
        print(f"## {title}\n  A = {la}\n  B = {lb}")
        for k, v in r.items():
            print(f"  {k}: {v}")
        print()
    return 0


# ------------------------------------------------------------------------ CLI

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("USAGE")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--mafft", type=Path, default=REPO / "mafft-upstream/scripts/mafft",
                        help="C mafft wrapper (default: in-tree build)")
        sp.add_argument("--mafft-binaries", type=Path,
                        default=REPO / "mafft-upstream/binaries",
                        help="MAFFT_BINARIES for the C wrapper")
        sp.add_argument("--mafft-rs", type=Path, default=REPO / "target/release/mafft-rs")
        sp.add_argument("--mafft-rs-nofma", type=Path,
                        default=REPO / "target/nofma/release/mafft-rs",
                        help="mafft-rs built with --features fp-contract-none")

    r = sub.add_parser("run", help="align (cached) and score")
    common(r)
    r.add_argument("--balibase", default=os.environ.get("BALIBASE_DIR"),
                   help="BAliBASE 3 root (default: $BALIBASE_DIR)")
    r.add_argument("--modes", nargs="+", default=list(MODES), choices=list(MODES))
    r.add_argument("--filter", default="", help="substring filter on the input path")
    r.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    r.add_argument("--timeout", type=int, default=3600, help="seconds per alignment")
    r.add_argument("--cache", type=Path, default=REPO / "target/balibase_scores")
    r.add_argument("--output", type=Path, default=None)
    r.add_argument("--overwrite", action="store_true",
                   help="do not merge with an existing output file")
    r.add_argument("--no-retry", action="store_true",
                   help="do not retry cached non-timeout failures")
    r.add_argument("--quiet", action="store_true")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("summarize", help="Markdown summary of a results TSV")
    s.add_argument("tsv")
    s.set_defaults(func=cmd_summarize)

    c = sub.add_parser("counterexamples", help="the 712/738 and tied-trace cases")
    common(c)
    c.set_defaults(func=cmd_counterexamples)

    x = sub.add_parser("crossarch",
                       help="host C vs the other FP-contraction policy on BAliBASE")
    common(x)
    x.add_argument("--balibase", default=os.environ.get("BALIBASE_DIR"))
    x.add_argument("--modes", nargs="+", default=["fftns2"], choices=list(MODES))
    x.add_argument("--filter", default="")
    x.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    x.add_argument("--timeout", type=int, default=3600)
    x.add_argument("--cache", type=Path, default=REPO / "target/balibase_scores")
    x.add_argument("--output", type=Path, default=None)
    x.add_argument("--no-retry", action="store_true")
    x.set_defaults(func=cmd_crossarch)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
