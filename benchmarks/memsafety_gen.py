#!/usr/bin/env python3
"""Generate the edge-case input sets for the memory-safety comparison (issue #27).

Writes deterministic files (fixed seed) under `target/memsafety/`, each
set in its own directory with a `manifest.tsv` holding one row per run:

    case <TAB> input <TAB> flags

`flags` may reference other generated files as `{dir}/<name>`; the
harness (`benchmarks/memsafety.py`) expands `{dir}` to the absolute
directory of the set. Three sets are written:

* `adversarial/` -- one case per read path or consistency check that a C
  implementation has to get right without bounds checks: very long
  header lines, very long sequences, many sequences, bytes outside
  ASCII, NUL bytes, empty records, CRLF and missing final newlines,
  duplicate names, gap-only sequences, `--add` with an "alignment"
  whose rows differ in length (with and without `--keeplength`),
  `--seed` files with ragged rows, and malformed `--treein` files
  (out-of-range, zero and negative indices, re-used nodes, Newick
  instead of MAFFT's pair format, too few lines, an overlong line).
* `sample_add/` -- small `--add` / `--addfragments` / `--keeplength`
  sets cut from `mafft-upstream/test/sample` (existing alignment = the
  first 30 rows of the C FFT-NS-2 alignment, added = the last 6
  sequences, ungapped or as fragments).
* `mutants/` -- mechanically mutated copies of `mafft-upstream/test/sample`
  (byte substitutions, truncations, stray `>` / `-` / NUL insertions,
  deleted newlines, duplicated header lines), run under FFT-NS-2 only.

USAGE
-----

    benchmarks/memsafety_gen.py                 # writes target/memsafety/
    benchmarks/memsafety_gen.py --out DIR --mutants 400
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SEED = 27

AA = "ACDEFGHIKLMNPQRSTVWY"
NT = "ACGT"

SAMPLE = REPO / "mafft-upstream" / "test" / "sample"
SAMPLE_FFTNS2 = REPO / "mafft-upstream" / "test" / "sample.fftns2"


def wrap(seq: str, width: int = 60) -> str:
    if width <= 0:
        return seq + "\n"
    return "".join(seq[i : i + width] + "\n" for i in range(0, len(seq), width))


def fasta(records: list[tuple[str, str]], width: int = 60) -> bytes:
    return "".join(f">{name}\n{wrap(seq, width)}" for name, seq in records).encode()


def rand_seq(rng: random.Random, alphabet: str, n: int) -> str:
    return "".join(rng.choice(alphabet) for _ in range(n))


def mutate(rng: random.Random, seq: str, alphabet: str, rate: float) -> str:
    out = []
    for c in seq:
        r = rng.random()
        if r < rate / 3:
            continue
        if r < 2 * rate / 3:
            out.append(rng.choice(alphabet))
        elif r < rate:
            out.append(c + rng.choice(alphabet))
        else:
            out.append(c)
    return "".join(out)


def family(rng: random.Random, alphabet: str, n: int, length: int, rate: float = 0.15) -> list[str]:
    root = rand_seq(rng, alphabet, length)
    return [mutate(rng, root, alphabet, rate) for _ in range(n)]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            recs.append((line[1:], ""))
        elif recs:
            recs[-1] = (recs[-1][0], recs[-1][1] + line.strip())
    return recs


# --------------------------------------------------------------- adversarial


def adversarial(out: Path) -> list[tuple[str, str, str]]:
    rng = random.Random(SEED)
    files: dict[str, bytes] = {}
    cases: list[tuple[str, str, str]] = []

    def add(case: str, name: str, data: bytes | None, *flag_sets: str) -> None:
        if data is not None:
            files[name] = data
        for flags in flag_sets or ("",):
            cases.append((case, name, flags))

    prot = family(rng, AA, 6, 120)
    dna = family(rng, NT, 6, 300)
    prot_recs = [(f"p{i + 1}", s) for i, s in enumerate(prot)]
    dna_recs = [(f"d{i + 1}", s) for i, s in enumerate(dna)]

    add("control_protein", "control_protein.fa", fasta(prot_recs), "", "--localpair --maxiterate 2")
    add("control_dna", "control_dna.fa", fasta(dna_recs), "", "--nuc")

    # Header lengths around C's fixed buffers (B = 256, 1000, N = 5000000
    # bytes) and far beyond them.
    for n in (255, 256, 257, 999, 1000, 1001, 10_000, 1_000_000):
        recs = [(("H" * n) if i == 0 else name, s) for i, (name, s) in enumerate(prot_recs)]
        add(f"long_header_{n}", f"long_header_{n}.fa", fasta(recs))
    recs = [(f"{'x' * 5000}_{i}", s) for i, (_, s) in enumerate(prot_recs)]
    add("long_header_all_5000", "long_header_all_5000.fa", fasta(recs), "", "--treeout")
    recs = [(name + " " + "\t" * 300, s) for name, s in prot_recs]
    add("header_trailing_tabs", "header_trailing_tabs.fa", fasta(recs))

    # One very long sequence among short ones; one >100k-residue line.
    long_dna = rand_seq(rng, NT, 120_000)
    recs = [("long120k", long_dna)] + [
        (f"frag{i}", mutate(rng, long_dna[o : o + 400], NT, 0.05))
        for i, o in enumerate((1000, 50_000, 110_000))
    ]
    add("long_seq_120k", "long_seq_120k.fa", fasta(recs), "--nuc")
    add("long_seq_120k_oneline", "long_seq_120k_oneline.fa", fasta(recs, width=0), "--nuc")
    recs = [("longprot", rand_seq(rng, AA, 30_000))] + prot_recs[:3]
    add("long_seq_prot_30k", "long_seq_prot_30k.fa", fasta(recs))

    # Thousands of short sequences.
    short = family(rng, AA, 3000, 30, rate=0.2)
    add("many_short_3000", "many_short_3000.fa", fasta([(f"s{i}", s) for i, s in enumerate(short)]))
    tiny = [rand_seq(rng, NT, rng.randint(1, 4)) for _ in range(2000)]
    add("many_tiny_2000", "many_tiny_2000.fa", fasta([(f"t{i}", s) for i, s in enumerate(tiny)]), "--nuc")

    # Bytes outside ASCII.
    body = fasta(prot_recs)
    add(
        "nonascii_header_utf8",
        "nonascii_header_utf8.fa",
        body.replace(b">p1", ">p1 \u00e9\u00e8\u4e2d\u6587 \U0001F9EC".encode()).replace(b">p3", b">p3 \xff\xfe\x80"),
    )
    seqs = bytearray(body)
    for _ in range(12):
        i = rng.randrange(len(seqs))
        if seqs[i] not in b">\n" and seqs[i - 1] != ord(">"):
            seqs[i] = rng.randrange(0x80, 0x100)
    add("nonascii_seq", "nonascii_seq.fa", bytes(seqs), "", "--anysymbol")
    add("nul_in_seq", "nul_in_seq.fa", body.replace(b"\n", b"\nAC\x00DE", 3), "", "--anysymbol")
    add("nul_in_header", "nul_in_header.fa", body.replace(b">p2", b">p2\x00hidden"))
    add("nul_only", "nul_only.fa", b"\x00" * 64)
    add("binary_garbage", "binary_garbage.fa", bytes(rng.randrange(256) for _ in range(4096)))
    add("utf16le", "utf16le.fa", body.decode().encode("utf-16"))

    # Record structure.
    def recs_to(*rs: tuple[str, str]) -> bytes:
        return fasta(list(rs))

    p = prot
    add("empty_record_middle", "empty_record_middle.fa", recs_to(("a", p[0]), ("b", ""), ("c", p[1]), ("d", p[2])))
    add("empty_record_last", "empty_record_last.fa", recs_to(("a", p[0]), ("b", p[1]), ("c", "")))
    add("empty_record_first", "empty_record_first.fa", recs_to(("a", ""), ("b", p[1]), ("c", p[2])))
    add("all_empty_records", "all_empty_records.fa", b">a\n>b\n>c\n")
    add("empty_file", "empty_file.fa", b"")
    add("only_newlines", "only_newlines.fa", b"\n\n\n")
    add("header_no_newline", "header_no_newline.fa", b">lonely")
    add("single_seq", "single_seq.fa", recs_to(("only", p[0])))
    add("two_seqs_one_empty", "two_seqs_one_empty.fa", recs_to(("a", p[0]), ("b", "")))
    add("no_header", "no_header.fa", wrap(p[0]).encode() + wrap(p[1]).encode())
    add("text_before_header", "text_before_header.fa", b"stray text\n" + body)
    add("empty_name", "empty_name.fa", recs_to(("", p[0]), ("", p[1]), ("", p[2])))
    add("gt_inside_seq", "gt_inside_seq.fa", body.replace(b"\n", b">\n", 2))
    add("double_gt", "double_gt.fa", body.replace(b">p", b">>p"))
    add("single_residue_seqs", "single_residue_seqs.fa", recs_to(("a", "M"), ("b", "K"), ("c", "W")))
    add("one_residue_vs_long", "one_residue_vs_long.fa", recs_to(("a", "M"), ("b", p[0]), ("c", p[1])))

    # Line endings.
    add("crlf", "crlf.fa", body.replace(b"\n", b"\r\n"))
    add("cr_only", "cr_only.fa", body.replace(b"\n", b"\r"))
    mixed = bytearray()
    for k, line in enumerate(body.split(b"\n")):
        mixed += line + (b"\r\n", b"\n", b"\r")[k % 3]
    add("mixed_eol", "mixed_eol.fa", bytes(mixed))
    add("no_trailing_newline", "no_trailing_newline.fa", body.rstrip(b"\n"))
    add("crlf_no_trailing_newline", "crlf_no_trailing_newline.fa", body.replace(b"\n", b"\r\n").rstrip(b"\r\n"))
    add("whitespace_lines", "whitespace_lines.fa", body.replace(b"\n>", b"\n \t \n\n>"))

    # Names and gaps.
    add("duplicate_names", "duplicate_names.fa", fasta([("same", s) for s in prot]))
    add(
        "duplicate_names_dna_adjust",
        "duplicate_names_dna.fa",
        fasta([("same", s) for s in dna]),
        "--adjustdirection --nuc",
    )
    add("gap_only_one", "gap_only_one.fa", recs_to(("a", p[0]), ("gaps", "-" * 120), ("c", p[1])))
    add("gap_only_all", "gap_only_all.fa", recs_to(("a", "-" * 50), ("b", "-" * 50), ("c", "---")))
    add("star_only_one", "star_only_one.fa", recs_to(("a", p[0]), ("stars", "*" * 40), ("c", p[1])))
    add("dot_only_dna", "dot_only_dna.fa", recs_to(("a", dna[0]), ("dots", "." * 40)), "--nuc", "--anysymbol --nuc")
    add("digits_only", "digits_only.fa", recs_to(("a", p[0]), ("digits", "1234567890" * 5)))

    # --add with mismatched lengths.
    existing = [(f"e{i}", s) for i, s in enumerate(family(rng, AA, 4, 100, rate=0.0))]
    ragged = [(n, s[: 100 - 7 * k] + "-" * (k % 2)) for k, (n, s) in enumerate(existing)]
    files["add_ragged_existing.fa"] = fasta(ragged)
    files["add_new.fa"] = fasta([("n1", mutate(rng, existing[0][1], AA, 0.1)), ("n2", existing[1][1][10:70])])
    files["add_new_long.fa"] = fasta([("nlong", existing[0][1] * 5)])
    files["add_new_empty.fa"] = b""
    files["add_new_gaponly.fa"] = fasta([("ngap", "-" * 30)])
    files["add_new_nonascii.fa"] = fasta([("nx", existing[2][1][:50])]).replace(b"A", b"\xc3\xa9")
    aligned = [(n, s) for n, s in existing]
    files["add_existing_ok.fa"] = fasta(aligned)
    for flags in ("--add {dir}/add_new.fa", "--add {dir}/add_new.fa --keeplength", "--addfragments {dir}/add_new.fa"):
        cases.append(("add_ragged_existing", "add_ragged_existing.fa", flags))
    for newfile in ("add_new_long", "add_new_empty", "add_new_gaponly", "add_new_nonascii"):
        for kl in ("", " --keeplength"):
            cases.append((f"{newfile}", "add_existing_ok.fa", f"--add {{dir}}/{newfile}.fa{kl}"))
    cases.append(("add_new_self", "add_existing_ok.fa", "--add {dir}/add_existing_ok.fa"))
    cases.append(("add_to_single", "single_seq.fa", "--add {dir}/add_new.fa --keeplength"))

    # --seed with ragged rows.
    seed = [(f"sd{i}", s) for i, s in enumerate(family(rng, AA, 3, 80, rate=0.0))]
    files["seed_ok.fa"] = fasta(seed)
    files["seed_ragged.fa"] = fasta([(n, s[: 80 - 9 * k]) for k, (n, s) in enumerate(seed)])
    files["seed_single.fa"] = fasta(seed[:1])
    files["seed_gaponly.fa"] = fasta(seed[:2] + [("sdgap", "-" * 80)])
    files["seed_empty.fa"] = b""
    for s in ("seed_ok", "seed_ragged", "seed_single", "seed_gaponly", "seed_empty"):
        cases.append((s, "control_protein.fa", f"--seed {{dir}}/{s}.fa"))
    cases.append(("seed_ragged_twice", "control_protein.fa", "--seed {dir}/seed_ragged.fa --seed {dir}/seed_ok.fa"))

    # --treein: MAFFT pair format, one merge per line, 1-based, i < j.
    trees = {
        "tree_ok": "1 2 0.1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
        "tree_zero_index": "0 2 0.1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
        "tree_negative_index": "-5 2 0.1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
        "tree_out_of_range": "1 99 0.1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
        "tree_reused_node": "1 2 0.1 0.1\n1 2 0.1 0.1\n1 2 0.1 0.1\n1 2 0.1 0.1\n1 2 0.1 0.1\n",
        "tree_merged_away": "1 2 0.1 0.1\n2 3 0.1 0.1\n2 4 0.1 0.1\n2 5 0.1 0.1\n2 6 0.1 0.1\n",
        "tree_equal_pair": "2 2 0.1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
        "tree_newick": "((p1:0.1,p2:0.1):0.2,((p3:0.1,p4:0.1):0.1,(p5:0.1,p6:0.1):0.1):0.2);\n",
        "tree_too_few_lines": "1 2 0.1 0.1\n3 4 0.1 0.1\n",
        "tree_empty": "",
        "tree_overlong_line": "1 2 " + "0" * 5000 + ".1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
        "tree_nan_lengths": "1 2 nan nan\n3 4 inf -inf\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 1e308 1e308\n",
        "tree_huge_index": "1 2147483647 0.1 0.1\n3 4 0.1 0.1\n1 3 0.2 0.2\n5 6 0.1 0.1\n1 5 0.3 0.3\n",
    }
    for name, text in trees.items():
        files[f"{name}.txt"] = text.encode()
        cases.append((name, "control_protein.fa", f"--treein {{dir}}/{name}.txt"))
    cases.append(("tree_zero_index_fftnsi", "control_protein.fa", "--treein {dir}/tree_zero_index.txt --maxiterate 2"))

    out.mkdir(parents=True, exist_ok=True)
    for name, data in files.items():
        (out / name).write_bytes(data)
    return cases


# ---------------------------------------------------------------- sample_add


def sample_add(out: Path) -> list[tuple[str, str, str]]:
    rng = random.Random(SEED + 1)
    aln = read_fasta(SAMPLE_FFTNS2)
    raw = read_fasta(SAMPLE)
    rows = aln[:30]
    keep = [j for j in range(len(rows[0][1])) if any(s[j] != "-" for _, s in rows)]
    existing = [(n, "".join(s[j] for j in keep)) for n, s in rows]
    new = [(n, s.replace("-", "")) for n, s in raw[30:]]
    frags = []
    for n, s in new:
        o = rng.randrange(0, max(1, len(s) - 120))
        frags.append((n + "_frag", s[o : o + rng.randint(40, 120)]))
    out.mkdir(parents=True, exist_ok=True)
    (out / "existing30.fa").write_bytes(fasta(existing))
    (out / "new6.fa").write_bytes(fasta(new))
    (out / "frag6.fa").write_bytes(fasta(frags))
    cases = []
    for flags in (
        "--add {dir}/new6.fa",
        "--add {dir}/new6.fa --keeplength",
        "--add {dir}/new6.fa --maxiterate 2",
        "--addfragments {dir}/frag6.fa",
        "--addfragments {dir}/frag6.fa --keeplength",
        "--addfragments {dir}/frag6.fa --localpair --maxiterate 2",
        "--add {dir}/new6.fa --reorder",
    ):
        cases.append(("sample_add", "existing30.fa", flags))
    return cases


# ------------------------------------------------------------------ mutants

POOL = [0, ord(">"), ord("-"), ord("\n"), ord("\r"), ord(" "), ord("*"), ord("."), ord("="), 0x80, 0xFF]


def mutants(out: Path, n: int) -> list[tuple[str, str, str]]:
    rng = random.Random(SEED + 2)
    base = SAMPLE.read_bytes()
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    ops = ("subst", "truncate", "insert_gt", "insert_gap", "insert_nul", "del_newline", "dup_header", "gt_linestart")
    for k in range(n):
        op = ops[k % len(ops)]
        b = bytearray(base)
        if op == "subst":
            for _ in range(rng.randint(1, 8)):
                b[rng.randrange(len(b))] = rng.choice(POOL) if rng.random() < 0.6 else rng.randrange(256)
        elif op == "truncate":
            b = b[: rng.randrange(len(b))]
        elif op in ("insert_gt", "insert_gap", "insert_nul"):
            ch = {"insert_gt": b">", "insert_gap": b"-", "insert_nul": b"\x00"}[op]
            for _ in range(rng.randint(1, 5)):
                i = rng.randrange(len(b) + 1)
                b[i:i] = ch * rng.randint(1, 3)
        elif op == "del_newline":
            nl = [i for i, c in enumerate(b) if c == 10]
            for i in sorted(rng.sample(nl, rng.randint(1, 6)), reverse=True):
                del b[i]
        elif op == "dup_header":
            lines = bytes(b).split(b"\n")
            hs = [i for i, line in enumerate(lines) if line.startswith(b">")]
            i = rng.choice(hs)
            lines.insert(i, lines[i])
            b = bytearray(b"\n".join(lines))
        elif op == "gt_linestart":
            lines = bytes(b).split(b"\n")
            for _ in range(rng.randint(1, 4)):
                i = rng.randrange(len(lines))
                lines[i] = b">" + lines[i]
            b = bytearray(b"\n".join(lines))
        name = f"mut{k:04d}_{op}.fa"
        (out / name).write_bytes(bytes(b))
        cases.append((name[:-3], name, ""))
    return cases


def write_manifest(out: Path, cases: list[tuple[str, str, str]]) -> None:
    with open(out / "manifest.tsv", "w") as f:
        f.write("case\tinput\tflags\n")
        for row in cases:
            f.write("\t".join(row) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=REPO / "target" / "memsafety")
    ap.add_argument("--mutants", type=int, default=400, help="number of mutated sample copies")
    a = ap.parse_args()
    for sub, fn in (
        ("adversarial", adversarial),
        ("sample_add", sample_add),
        ("mutants", lambda d: mutants(d, a.mutants)),
    ):
        d = a.out / sub
        cases = fn(d)
        write_manifest(d, cases)
        print(f"{d}: {len(cases)} runs")


if __name__ == "__main__":
    main()
