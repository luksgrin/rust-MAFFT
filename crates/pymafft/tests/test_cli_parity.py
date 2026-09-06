"""Option-by-option parity of pymafft against the `mafft-rs` command line.

Every keyword of `pymafft.align` maps to one CLI flag and is applied by the
CLI's own flag layer (`mafft_rs::run_from_seqs`), so for each option set the
Python result must be byte-identical to what the binary prints on stdout for
the same fixture and flags. The same holds for errors (`MafftError.code` /
`.message` are the exit code and stderr text) and for progress lines.

The binary is the release build of THIS checkout (`cargo build --release -p
mafft-rs`); the module skips if it is missing.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

import pymafft


REPO_ROOT = Path(__file__).resolve().parents[3]
BIN_FIXTURES = REPO_ROOT / "crates" / "mafft-bin" / "tests" / "fixtures"
CORE_FIXTURES = REPO_ROOT / "crates" / "mafft-core" / "tests" / "fixtures"
UPSTREAM_TEST = REPO_ROOT / "mafft-upstream" / "test"

DNA_PAIR = BIN_FIXTURES / "dna_pair_gapscale_min.fa"
REFINE = BIN_FIXTURES / "refine_njob2_min.fa"
MTB = BIN_FIXTURES / "mtb_cds_120x1400.fa"
ILLEGAL_PROTEIN = BIN_FIXTURES / "input_handling" / "prot_unusual_UOJBZX.fa"
SAMPLE = UPSTREAM_TEST / "sample"
SAMPLE_FFTNS2_EXPECTED = UPSTREAM_TEST / "sample.fftns2"
SAMPLE_FIRST14 = CORE_FIXTURES / "sample.first14.fa"

PANAROO = dict(strategy="auto", adjust_direction=True, threads=1, seq_type="nuc")
PANAROO_FLAGS = ["--auto", "--adjustdirection", "--thread", "1", "--nuc"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _find_mafft_rs() -> str | None:
    local = REPO_ROOT / "target" / "release" / "mafft-rs"
    if local.exists():
        return str(local)
    return shutil.which("mafft-rs")


@pytest.fixture(scope="module")
def mafft_rs() -> str:
    binary = _find_mafft_rs()
    if binary is None:
        pytest.skip("mafft-rs binary not found (run `cargo build --release -p mafft-rs`)")
    return binary


def _cli(binary: str, flags: list[str], path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [binary, *flags, str(path)], capture_output=True, check=False, timeout=600
    )


def _unwrap_fasta(text: str) -> str:
    """Join the CLI's 60-column FASTA into `to_fasta()`'s one-line-per-row form."""
    records: list[str] = []
    seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if seq:
                records.append("".join(seq))
                seq = []
            records.append(line)
        else:
            seq.append(line)
    if seq:
        records.append("".join(seq))
    return "\n".join(records)


def _cli_fasta(binary: str, flags: list[str], path: Path) -> str:
    proc = _cli(binary, flags, path)
    assert proc.returncode == 0, proc.stderr.decode()
    return _unwrap_fasta(proc.stdout.decode())


def _read_records(path: Path) -> list[tuple[str, str]]:
    """Raw FASTA records — residues verbatim, so the in-memory path has to
    normalise them exactly as the reader does for the CLI."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].rstrip()
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _needs(*paths: Path):
    missing = [p.name for p in paths if not p.exists()]
    return pytest.mark.skipif(bool(missing), reason=f"fixture missing: {missing}")


# ---------------------------------------------------------------------------
# Keyword -> flag mapping
# ---------------------------------------------------------------------------


class TestFlagMapping:
    """The option table itself, before any alignment runs."""

    def _flags(self, **kw):
        base = dict(
            seq_type=None, adjust_direction=False, threads=0, reorder=False,
            retree=None, scoring=None, gap_open=None, gap_extend=None,
            quiet=True, progress=None,
        )
        base.update(kw)
        strategy = base.pop("strategy", "fftns2")
        maxiterate = base.pop("maxiterate", 0)
        return pymafft._flags(strategy, maxiterate, **base)

    def test_defaults_are_fftns2_quiet(self):
        assert self._flags() == ["--quiet"]

    def test_panaroo_set(self):
        assert self._flags(**PANAROO) == [*PANAROO_FLAGS, "--quiet"]

    def test_strategy_defaults_match_previous_pymafft(self):
        assert self._flags(strategy="fftnsi") == ["--maxiterate", "100", "--quiet"]
        assert self._flags(strategy="ginsi") == ["--globalpair", "--maxiterate", "1000", "--quiet"]
        assert self._flags(strategy="linsi") == ["--localpair", "--maxiterate", "1000", "--quiet"]
        assert self._flags(strategy="einsi") == ["--genafpair", "--maxiterate", "1000", "--quiet"]
        assert self._flags(strategy="linsi", maxiterate=2) == ["--localpair", "--maxiterate", "2", "--quiet"]
        # fftns2 has always ignored maxiterate.
        assert self._flags(strategy="fftns2", maxiterate=50) == ["--quiet"]
        assert self._flags(strategy="auto", maxiterate=50) == ["--auto", "--maxiterate", "50", "--quiet"]

    def test_every_keyword(self):
        flags = self._flags(
            strategy="auto", seq_type="amino", adjust_direction="accurately",
            threads=4, reorder=True, retree=1, scoring="bl45", gap_open=2.0,
            gap_extend=0.5,
        )
        assert flags == [
            "--auto", "--adjustdirectionaccurately", "--thread", "4", "--amino",
            "--reorder", "--retree", "1", "--bl", "45", "--op", "2.0", "--ep", "0.5",
            "--quiet",
        ]

    def test_scoring_spellings(self):
        assert self._flags(scoring="BL62")[:2] == ["--bl", "62"]
        assert self._flags(scoring="blosum80")[:2] == ["--bl", "80"]
        assert self._flags(scoring="jtt200")[:2] == ["--jtt", "200"]
        assert self._flags(scoring="tm100")[:2] == ["--tm", "100"]

    def test_progress_disables_quiet(self):
        assert self._flags(progress=print) == []
        assert self._flags(quiet=False) == []

    @pytest.mark.parametrize(
        "kw, exc, match",
        [
            (dict(strategy="bogus"), ValueError, "unknown strategy"),
            (dict(seq_type="dna"), ValueError, "seq_type"),
            (dict(adjust_direction="yes"), ValueError, "adjust_direction"),
            (dict(scoring="pam250"), ValueError, "scoring"),
            (dict(threads=-1), ValueError, "threads"),
            (dict(threads=1.5), TypeError, "threads"),
            (dict(retree=-1), ValueError, "retree"),
            (dict(maxiterate=-1), ValueError, "maxiterate"),
        ],
    )
    def test_invalid_values(self, kw, exc, match):
        with pytest.raises(exc, match=match):
            self._flags(**kw)


# ---------------------------------------------------------------------------
# Alignment parity, option by option
# ---------------------------------------------------------------------------

# (fixture, pymafft keywords, equivalent CLI flags). The CLI always gets
# `--quiet` because pymafft's default is quiet=True; it does not affect the
# alignment.
CASES = [
    pytest.param(DNA_PAIR, {}, [], id="dna_pair-default"),
    pytest.param(DNA_PAIR, dict(strategy="auto"), ["--auto"], id="dna_pair-auto"),
    pytest.param(DNA_PAIR, dict(strategy="linsi", maxiterate=1000),
                 ["--localpair", "--maxiterate", "1000"], id="dna_pair-linsi"),
    pytest.param(DNA_PAIR, dict(strategy="einsi", maxiterate=2),
                 ["--genafpair", "--maxiterate", "2"], id="dna_pair-einsi"),
    pytest.param(DNA_PAIR, dict(seq_type="nuc"), ["--nuc"], id="dna_pair-nuc"),
    pytest.param(DNA_PAIR, dict(adjust_direction=True), ["--adjustdirection"],
                 id="dna_pair-adjustdirection"),
    pytest.param(DNA_PAIR, dict(adjust_direction="accurately"),
                 ["--adjustdirectionaccurately"], id="dna_pair-adjustdirectionaccurately"),
    pytest.param(DNA_PAIR, dict(gap_open=3.0, gap_extend=0.5), ["--op", "3.0", "--ep", "0.5"],
                 id="dna_pair-op-ep"),
    pytest.param(DNA_PAIR, dict(retree=1), ["--retree", "1"], id="dna_pair-retree1"),
    pytest.param(REFINE, dict(strategy="fftnsi", maxiterate=2), ["--maxiterate", "2"],
                 id="refine-fftnsi2"),
    pytest.param(REFINE, dict(strategy="auto", threads=1), ["--auto", "--thread", "1"],
                 id="refine-auto-thread1"),
    pytest.param(REFINE, dict(strategy="ginsi", maxiterate=2, threads=1),
                 ["--globalpair", "--maxiterate", "2", "--thread", "1"], id="refine-ginsi-thread1"),
    pytest.param(REFINE, dict(reorder=True), ["--reorder"], id="refine-reorder"),
    pytest.param(MTB, PANAROO, PANAROO_FLAGS, id="mtb-panaroo"),
    pytest.param(MTB, dict(**PANAROO, reorder=True), [*PANAROO_FLAGS, "--reorder"],
                 id="mtb-panaroo-reorder"),
    pytest.param(SAMPLE, dict(scoring="bl50", retree=2), ["--bl", "50", "--retree", "2"],
                 id="sample-bl50"),
    pytest.param(SAMPLE, dict(scoring="jtt200", strategy="fftnsi", maxiterate=2),
                 ["--jtt", "200", "--maxiterate", "2"], id="sample-jtt200-fftnsi2"),
    pytest.param(SAMPLE, dict(scoring="tm100"), ["--tm", "100"], id="sample-tm100"),
    pytest.param(SAMPLE, dict(seq_type="amino"), ["--amino"], id="sample-amino"),
    pytest.param(SAMPLE, dict(strategy="auto"), ["--auto"], id="sample-auto"),
    pytest.param(SAMPLE, dict(reorder=True), ["--reorder"], id="sample-reorder"),
    pytest.param(SAMPLE, dict(gap_open=1.0, gap_extend=0.2), ["--op", "1.0", "--ep", "0.2"],
                 id="sample-op-ep"),
    pytest.param(SAMPLE_FIRST14, dict(strategy="linsi", maxiterate=2, scoring="bl45"),
                 ["--localpair", "--maxiterate", "2", "--bl", "45"], id="first14-linsi-bl45"),
]


class TestOptionParity:
    """`align_file` and `align` (in-memory records) both equal the CLI."""

    @pytest.mark.parametrize("fixture, kwargs, flags", CASES)
    def test_align_file_matches_cli(self, mafft_rs, fixture, kwargs, flags):
        if not fixture.exists():
            pytest.skip(f"fixture missing: {fixture.name}")
        expected = _cli_fasta(mafft_rs, [*flags, "--quiet"], fixture)
        result = pymafft.align_file(str(fixture), **kwargs)
        assert result.to_fasta() == expected

    @pytest.mark.parametrize("fixture, kwargs, flags", CASES)
    def test_align_records_matches_cli(self, mafft_rs, fixture, kwargs, flags):
        if not fixture.exists():
            pytest.skip(f"fixture missing: {fixture.name}")
        expected = _cli_fasta(mafft_rs, [*flags, "--quiet"], fixture)
        result = pymafft.align(_read_records(fixture), **kwargs)
        assert result.to_fasta() == expected

    @_needs(MTB)
    def test_panaroo_set_over_plain_strings(self, mafft_rs):
        """The panaroo call shape: bare uppercase strings, no names."""
        expected = _cli_fasta(mafft_rs, [*PANAROO_FLAGS, "--quiet"], MTB)
        records = _read_records(MTB)
        result = pymafft.align([seq.upper() for _, seq in records], **PANAROO)
        # Names are auto-generated; the rows must still be identical.
        assert [s for _, s in result.to_tuples()] == [
            s for _, s in pymafft.align(records, **PANAROO).to_tuples()
        ]
        assert [s for _, s in result.to_tuples()] == [
            line for line in expected.splitlines() if not line.startswith(">")
        ]

    @_needs(SAMPLE)
    def test_fasta_string_matches_cli(self, mafft_rs):
        expected = _cli_fasta(mafft_rs, ["--bl", "50", "--retree", "2", "--quiet"], SAMPLE)
        result = pymafft.align_fasta_string(SAMPLE.read_text(), scoring="bl50", retree=2)
        assert result.to_fasta() == expected


# ---------------------------------------------------------------------------
# `run()` escape hatch
# ---------------------------------------------------------------------------


class TestRun:
    @_needs(MTB)
    def test_run_equals_align_for_equivalent_flags(self):
        records = _read_records(MTB)
        via_run = pymafft.run([*PANAROO_FLAGS, "--quiet"], records)
        via_align = pymafft.align(records, **PANAROO)
        assert via_run.to_fasta() == via_align.to_fasta()

    @_needs(SAMPLE)
    def test_run_matches_cli_for_flag_without_keyword(self, mafft_rs):
        flags = ["--nofft", "--maxiterate", "2", "--quiet"]
        expected = _cli_fasta(mafft_rs, flags, SAMPLE)
        assert pymafft.run(flags, _read_records(SAMPLE)).to_fasta() == expected

    def test_run_rejects_bare_string(self):
        with pytest.raises(TypeError):
            pymafft.run("--auto", ["ACGT", "ACGA"])

    def test_run_rejects_input_path(self):
        with pytest.raises(pymafft.MafftError) as info:
            pymafft.run(["--quiet", "in.fa"], ["ACGT", "ACGA"])
        assert info.value.code == 1
        assert "INPUT" in info.value.message

    def test_unknown_flag_matches_cli(self, mafft_rs):
        proc = subprocess.run([mafft_rs, "--no-such-flag"], capture_output=True, check=False)
        with pytest.raises(pymafft.MafftError) as info:
            pymafft.run(["--no-such-flag"], ["ACGT", "ACGA"])
        assert info.value.code == proc.returncode
        assert info.value.message.strip() == proc.stderr.decode().strip()


# ---------------------------------------------------------------------------
# Error parity
# ---------------------------------------------------------------------------


class TestErrorParity:
    @_needs(ILLEGAL_PROTEIN)
    def test_illegal_residue_quiet(self, mafft_rs):
        proc = _cli(mafft_rs, ["--quiet"], ILLEGAL_PROTEIN)
        assert proc.returncode == 1
        with pytest.raises(pymafft.MafftError) as info:
            pymafft.align_file(str(ILLEGAL_PROTEIN))
        err = info.value
        assert isinstance(err, ValueError)
        assert err.code == proc.returncode
        assert err.message + "\n" == proc.stderr.decode()
        assert str(err) == err.message == "Illegal character U"

    @_needs(ILLEGAL_PROTEIN)
    def test_illegal_residue_with_banner(self, mafft_rs):
        proc = _cli(mafft_rs, [], ILLEGAL_PROTEIN)
        assert proc.returncode == 1
        with pytest.raises(pymafft.MafftError) as info:
            pymafft.align_file(str(ILLEGAL_PROTEIN), quiet=False)
        assert info.value.code == 1
        assert info.value.message + "\n" == proc.stderr.decode()
        assert "Alphabet 'U' is unknown" in info.value.message

    @_needs(ILLEGAL_PROTEIN)
    def test_illegal_residue_in_memory(self, mafft_rs):
        proc = _cli(mafft_rs, ["--quiet"], ILLEGAL_PROTEIN)
        with pytest.raises(pymafft.MafftError) as info:
            pymafft.align(_read_records(ILLEGAL_PROTEIN))
        assert (info.value.code, info.value.message + "\n") == (1, proc.stderr.decode())

    def test_nucleotide_dot_is_illegal(self, mafft_rs, tmp_path):
        fa = tmp_path / "dot.fa"
        fa.write_text(">a\nacgt.acgt\n>b\nacgtacgt\n")
        proc = _cli(mafft_rs, ["--quiet"], fa)
        assert proc.returncode == 1
        with pytest.raises(ValueError) as info:  # MafftError is a ValueError
            pymafft.align([("a", "acgt.acgt"), ("b", "acgtacgt")])
        assert isinstance(info.value, pymafft.MafftError)
        assert info.value.message + "\n" == proc.stderr.decode()

    def test_plain_value_errors_are_unchanged(self):
        with pytest.raises(ValueError, match="no sequences"):
            pymafft.align([])
        with pytest.raises(ValueError, match="error reading"):
            pymafft.align_file("/nonexistent/file.fa")


# ---------------------------------------------------------------------------
# Progress
# ---------------------------------------------------------------------------


class TestProgress:
    @_needs(SAMPLE)
    def test_callable_receives_cli_stderr_lines(self, mafft_rs):
        proc = _cli(mafft_rs, [], SAMPLE)
        lines: list[str] = []
        pymafft.align_file(str(SAMPLE), progress=lines.append)
        assert lines == proc.stderr.decode().splitlines()
        assert lines[0].startswith("mafft-rs v")

    @_needs(MTB)
    def test_callable_receives_cli_stderr_lines_panaroo(self, mafft_rs):
        proc = _cli(mafft_rs, PANAROO_FLAGS, MTB)
        lines: list[str] = []
        pymafft.align(_read_records(MTB), progress=lines.append, **PANAROO)
        assert lines == proc.stderr.decode().splitlines()

    def test_callable_with_run(self):
        lines: list[str] = []
        pymafft.run(["--auto"], ["ACGTACGTAC", "ACGTACGAC"], progress=lines.append)
        assert any("strategy" in line for line in lines)

    def test_quiet_default_sends_nothing_to_stderr(self, capfd):
        pymafft.align(["ACDEFGHIK", "ACDEFHIK"])
        assert capfd.readouterr().err == ""

    def test_quiet_false_goes_to_stderr(self, capfd):
        pymafft.align(["ACDEFGHIK", "ACDEFHIK"], quiet=False)
        err = capfd.readouterr().err
        assert "strategy: FFT-NS-2" in err

    def test_callable_silences_stderr(self, capfd):
        seen: list[str] = []
        pymafft.align(["ACDEFGHIK", "ACDEFHIK"], progress=seen.append)
        assert capfd.readouterr().err == ""
        assert seen

    def test_callback_exception_propagates(self):
        def boom(_msg: str) -> None:
            raise RuntimeError("stop")

        with pytest.raises(RuntimeError, match="stop"):
            pymafft.align(["ACDEFGHIK", "ACDEFHIK"], progress=boom)


# ---------------------------------------------------------------------------
# Regression: the previous defaults still give the same bytes
# ---------------------------------------------------------------------------


class TestDefaultsRegression:
    @_needs(SAMPLE, SAMPLE_FFTNS2_EXPECTED)
    def test_default_align_file_equals_c_fixture(self):
        expected = _unwrap_fasta(SAMPLE_FFTNS2_EXPECTED.read_text())
        assert pymafft.align_file(str(SAMPLE)).to_fasta() == expected

    @_needs(SAMPLE, SAMPLE_FFTNS2_EXPECTED)
    def test_default_align_records_equals_align_file(self):
        from_file = pymafft.align_file(str(SAMPLE))
        from_records = pymafft.align(_read_records(SAMPLE))
        from_string = pymafft.align_fasta_string(SAMPLE.read_text())
        assert from_records.to_fasta() == from_file.to_fasta() == from_string.to_fasta()

    @_needs(SAMPLE)
    def test_default_equals_cli_default(self, mafft_rs):
        assert pymafft.align_file(str(SAMPLE)).to_fasta() == _cli_fasta(mafft_rs, ["--quiet"], SAMPLE)
