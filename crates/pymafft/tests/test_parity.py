"""Parity tests for pymafft.

These tests cross-validate that the Python bindings produce the SAME alignment
output as:
  1. The `mafft-rs` binary (proves the Python wrapper doesn't change anything).
  2. The system `mafft` binary (proves byte-identity with C MAFFT 7.526
     extends to the Python path).

The tests are skipped if `mafft-rs` or `mafft` aren't on PATH, so the suite
stays green on machines that only have one or the other.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

import pymafft


REPO_ROOT = Path(__file__).resolve().parents[3]
UPSTREAM_DIR = REPO_ROOT / "mafft-upstream" / "test"
FIXTURES_DIR = REPO_ROOT / "crates" / "mafft-core" / "tests" / "fixtures"

SAMPLE = UPSTREAM_DIR / "sample"
SAMPLE_FIRST14 = FIXTURES_DIR / "sample.first14.fa"
SAMPLE_FIRST15 = FIXTURES_DIR / "sample.first15.fa"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _find_mafft_rs() -> str | None:
    """Locate the `mafft-rs` binary. Prefer the local release build."""
    local = REPO_ROOT / "target" / "release" / "mafft-rs"
    if local.exists():
        return str(local)
    return shutil.which("mafft-rs")


def _find_c_mafft() -> str | None:
    """Locate the system `mafft` binary (C MAFFT)."""
    return shutil.which("mafft")


def _run_binary(binary: str, args: list[str], stdin: bytes | None = None) -> str:
    """Run a binary, return stdout. Raises on non-zero exit.

    C MAFFT's shell wrapper requires all flags BEFORE the input path:
        mafft --maxiterate 100 input.fa     # works
        mafft input.fa --maxiterate 100     # error
    Tests should pass `args` already ordered as [flags..., input_path].
    """
    result = subprocess.run(
        [binary, *args],
        input=stdin,
        capture_output=True,
        check=True,
    )
    return result.stdout.decode("utf-8")


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA into ordered list of (name, sequence) tuples."""
    entries: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                entries.append((name, "".join(chunks)))
            name = line[1:].rstrip()
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        entries.append((name, "".join(chunks)))
    return entries


def _normalize_alignment(entries: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """Uppercase residues, strip whitespace. Mirrors how pymafft processes its input."""
    return [(n, s.upper()) for n, s in entries]


# ---------------------------------------------------------------------------
# Round-trip tests (no external binary needed)
# ---------------------------------------------------------------------------


class TestRoundTrip:
    """`align` and `align_fasta_string` of the same data give the same alignment."""

    def test_strings_vs_tuples_same_alignment(self):
        seqs = ["ACDEFGHIK", "ACDEFHIK", "ACDHIK"]
        a = pymafft.align(seqs)
        b = pymafft.align([(f"seq_{i+1}", s) for i, s in enumerate(seqs)])
        assert a.to_tuples() == b.to_tuples()

    def test_fasta_string_vs_tuples(self):
        seqs = [("alpha", "ACDEFGHIK"), ("beta", "ACDEFHIK"), ("gamma", "ACDHIK")]
        from_tuples = pymafft.align(seqs)
        fasta_in = "\n".join(f">{n}\n{s}" for n, s in seqs)
        from_string = pymafft.align_fasta_string(fasta_in)
        assert from_tuples.to_tuples() == from_string.to_tuples()

    def test_to_fasta_reparses_to_same_alignment(self):
        result = pymafft.align(
            [("a", "ACDEFGHIK"), ("b", "ACDEFHIK"), ("c", "ACDHIK")]
        )
        fasta = result.to_fasta()
        parsed = _parse_fasta(fasta)
        assert parsed == result.to_tuples()

    def test_deterministic_across_calls(self):
        seqs = ["MNGTEGDNFY", "MNGTEGPNFY", "MNGTEGINFY"]
        a = pymafft.align(seqs)
        b = pymafft.align(seqs)
        c = pymafft.align(seqs)
        assert a.to_tuples() == b.to_tuples() == c.to_tuples()


# ---------------------------------------------------------------------------
# Parity vs mafft-rs binary
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def mafft_rs():
    binary = _find_mafft_rs()
    if binary is None:
        pytest.skip("mafft-rs binary not found (run `cargo build --release`)")
    return binary


@pytest.fixture(scope="module")
def c_mafft():
    binary = _find_c_mafft()
    if binary is None:
        pytest.skip("system mafft (C) binary not found")
    return binary


class TestParityWithMafftRs:
    """pymafft and the mafft-rs CLI produce the same alignment."""

    @pytest.mark.skipif(not SAMPLE.exists(), reason="upstream submodule missing")
    def test_default_fftns2_36seq(self, mafft_rs):
        py_result = pymafft.align_file(str(SAMPLE))
        bin_text = _run_binary(mafft_rs, ["--quiet", str(SAMPLE)])
        bin_entries = _normalize_alignment(_parse_fasta(bin_text))
        assert py_result.to_tuples() == bin_entries

    @pytest.mark.skipif(not SAMPLE_FIRST14.exists(), reason="fixture missing")
    def test_fftnsi_first14(self, mafft_rs):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST14), strategy="fftnsi", maxiterate=100
        )
        bin_text = _run_binary(
            mafft_rs,
            ["--maxiterate", "100", "--quiet", str(SAMPLE_FIRST14)],
        )
        bin_entries = _normalize_alignment(_parse_fasta(bin_text))
        assert py_result.to_tuples() == bin_entries

    @pytest.mark.skipif(not SAMPLE_FIRST14.exists(), reason="fixture missing")
    def test_linsi_first14(self, mafft_rs):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST14), strategy="linsi", maxiterate=1000
        )
        bin_text = _run_binary(
            mafft_rs,
            [
                "--localpair", "--maxiterate", "1000", "--quiet",
                str(SAMPLE_FIRST14),
            ],
        )
        bin_entries = _normalize_alignment(_parse_fasta(bin_text))
        assert py_result.to_tuples() == bin_entries

    @pytest.mark.skipif(not SAMPLE_FIRST14.exists(), reason="fixture missing")
    def test_ginsi_first14(self, mafft_rs):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST14), strategy="ginsi", maxiterate=1000
        )
        bin_text = _run_binary(
            mafft_rs,
            [
                "--globalpair", "--maxiterate", "1000", "--quiet",
                str(SAMPLE_FIRST14),
            ],
        )
        bin_entries = _normalize_alignment(_parse_fasta(bin_text))
        assert py_result.to_tuples() == bin_entries

    @pytest.mark.skipif(not SAMPLE_FIRST14.exists(), reason="fixture missing")
    def test_einsi_first14(self, mafft_rs):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST14), strategy="einsi", maxiterate=1000
        )
        bin_text = _run_binary(
            mafft_rs,
            [
                "--genafpair", "--maxiterate", "1000", "--quiet",
                str(SAMPLE_FIRST14),
            ],
        )
        bin_entries = _normalize_alignment(_parse_fasta(bin_text))
        assert py_result.to_tuples() == bin_entries


# ---------------------------------------------------------------------------
# Parity vs C MAFFT
# ---------------------------------------------------------------------------


class TestParityWithCMafft:
    """pymafft and C MAFFT 7.526 produce the same alignment (mod case)."""

    @pytest.mark.skipif(not SAMPLE.exists(), reason="upstream submodule missing")
    def test_default_fftns2_36seq(self, c_mafft):
        py_result = pymafft.align_file(str(SAMPLE))
        c_text = _run_binary(c_mafft, ["--quiet", str(SAMPLE)])
        c_entries = _normalize_alignment(_parse_fasta(c_text))
        # Sequence content must match byte-for-byte; widths must match.
        py_tuples = py_result.to_tuples()
        assert len(py_tuples) == len(c_entries)
        assert all(len(p[1]) == len(c[1]) for p, c in zip(py_tuples, c_entries))
        for (pn, ps), (cn, cs) in zip(py_tuples, c_entries):
            assert pn == cn, f"name mismatch: {pn!r} vs {cn!r}"
            assert ps == cs, (
                f"sequence mismatch for {pn!r}:\n"
                f"  pymafft: {ps[:80]}...\n  C MAFFT: {cs[:80]}..."
            )

    @pytest.mark.skipif(not SAMPLE_FIRST15.exists(), reason="fixture missing")
    def test_fftnsi_first15(self, c_mafft):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST15), strategy="fftnsi", maxiterate=100
        )
        c_text = _run_binary(
            c_mafft,
            ["--maxiterate", "100", "--quiet", str(SAMPLE_FIRST15)],
        )
        c_entries = _normalize_alignment(_parse_fasta(c_text))
        py_tuples = py_result.to_tuples()
        for (pn, ps), (cn, cs) in zip(py_tuples, c_entries):
            assert ps == cs, f"divergence at {pn!r}"

    @pytest.mark.skipif(not SAMPLE_FIRST14.exists(), reason="fixture missing")
    def test_linsi_first14(self, c_mafft):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST14), strategy="linsi", maxiterate=1000
        )
        c_text = _run_binary(
            c_mafft,
            [
                "--localpair", "--maxiterate", "1000", "--quiet",
                str(SAMPLE_FIRST14),
            ],
        )
        c_entries = _normalize_alignment(_parse_fasta(c_text))
        py_tuples = py_result.to_tuples()
        for (pn, ps), (cn, cs) in zip(py_tuples, c_entries):
            assert ps == cs, f"divergence at {pn!r}"

    @pytest.mark.skipif(not SAMPLE_FIRST14.exists(), reason="fixture missing")
    def test_ginsi_first14(self, c_mafft):
        py_result = pymafft.align_file(
            str(SAMPLE_FIRST14), strategy="ginsi", maxiterate=1000
        )
        c_text = _run_binary(
            c_mafft,
            [
                "--globalpair", "--maxiterate", "1000", "--quiet",
                str(SAMPLE_FIRST14),
            ],
        )
        c_entries = _normalize_alignment(_parse_fasta(c_text))
        py_tuples = py_result.to_tuples()
        for (pn, ps), (cn, cs) in zip(py_tuples, c_entries):
            assert ps == cs, f"divergence at {pn!r}"


# ---------------------------------------------------------------------------
# Output-shape invariants
# ---------------------------------------------------------------------------


class TestOutputShape:
    """Invariants every alignment output must satisfy."""

    @pytest.fixture(scope="class")
    def result(self):
        return pymafft.align(
            [
                ("a", "MNGTEGDNFYVPFSNKTGVVRSPFEQPQYYLAEPWQF"),
                ("b", "MNGTEGPNFYVPFSNITGVVRSPFEQPQYYLAEPWQF"),
                ("c", "MNGTEGINFYVPMSNKTGVVRSPFEYPQYYLAEPWKY"),
                ("d", "MNGTEGKNFYVPMSNRTGLVRSPFEYPQYYLAEPWQF"),
            ]
        )

    def test_all_sequences_same_width(self, result):
        widths = {len(s) for s in result.sequences}
        assert len(widths) == 1, f"sequences have different widths: {widths}"

    def test_width_equals_first_sequence_length(self, result):
        assert result.width == len(result.sequences[0])

    def test_ungapped_equals_original_input(self, result):
        originals = [
            "MNGTEGDNFYVPFSNKTGVVRSPFEQPQYYLAEPWQF",
            "MNGTEGPNFYVPFSNITGVVRSPFEQPQYYLAEPWQF",
            "MNGTEGINFYVPMSNKTGVVRSPFEYPQYYLAEPWKY",
            "MNGTEGKNFYVPMSNRTGLVRSPFEYPQYYLAEPWQF",
        ]
        for original, seq in zip(originals, result.sequences):
            assert seq.ungapped() == original

    def test_indexable_and_iterable_agree(self, result):
        via_index = [result[i] for i in range(result.nseq)]
        via_iter = list(result)
        assert [s.name for s in via_index] == [s.name for s in via_iter]
        assert [s.sequence for s in via_index] == [s.sequence for s in via_iter]

    def test_fasta_roundtrips_via_align_fasta_string(self, result):
        fasta = result.to_fasta()
        reparsed = pymafft.align_fasta_string(fasta)
        # Aligning an already-aligned FASTA might collapse columns; but the
        # ungapped residues must be preserved.
        for original, reparsed_seq in zip(result.sequences, reparsed.sequences):
            assert original.ungapped() == reparsed_seq.ungapped()


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_single_sequence(self):
        result = pymafft.align(["ACDEFGHIK"])
        assert result.nseq == 1
        assert result.sequences[0].ungapped() == "ACDEFGHIK"

    def test_two_identical_sequences(self):
        result = pymafft.align(["ACDEFGHIK", "ACDEFGHIK"])
        assert result.sequences[0].sequence == result.sequences[1].sequence

    def test_very_short_sequences(self):
        result = pymafft.align(["AC", "AG", "AT"])
        assert result.nseq == 3
        for s in result.sequences:
            assert len(s.ungapped()) == 2

    def test_mixed_lengths(self):
        seqs = ["A" * 50, "A" * 10, "A" * 100]
        result = pymafft.align(seqs)
        assert result.width >= 100  # at least as wide as longest input
        for original, aligned in zip(seqs, result.sequences):
            assert aligned.ungapped() == original

    def test_case_passes_through_verbatim(self):
        # pymafft passes residues through to the engine as-is (no auto-
        # uppercase). The scoring matrix only has entries for uppercase
        # protein residues — lowercase chars are scored as unknowns, so
        # lowercase input may align differently from uppercase input.
        # Just verify case is preserved end-to-end (`align` doesn't mutate
        # the case of residues that survive alignment).
        upper = pymafft.align(["ACDEFGHIK", "ACDEFHIK"])
        assert upper.sequences[0].ungapped() == "ACDEFGHIK"

        lower = pymafft.align(["acdefghik", "acdefhik"])
        assert lower.sequences[0].ungapped() == "acdefghik"
