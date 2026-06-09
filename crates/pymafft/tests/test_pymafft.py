"""Test suite for pymafft Python bindings."""

import os
import pytest
import pymafft


# Path to MAFFT upstream test data
TEST_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "mafft-upstream", "test"
)
SAMPLE_PATH = os.path.join(TEST_DIR, "sample")
SAMPLERNA_PATH = os.path.join(TEST_DIR, "samplerna")


# ---------------------------------------------------------------------------
# Basic API tests
# ---------------------------------------------------------------------------


class TestAlign:
    """Test the align() function."""

    def test_list_of_strings(self):
        result = pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])
        assert result.nseq == 3
        assert result.width > 0
        # All sequences same width
        for seq in result.sequences:
            assert len(seq) == result.width

    def test_list_of_tuples(self):
        seqs = [("s1", "ACDEFGHIK"), ("s2", "ACDEFHIK"), ("s3", "ACDHIK")]
        result = pymafft.align(seqs)
        assert result.nseq == 3
        assert result.sequences[0].name == "s1"
        assert result.sequences[1].name == "s2"

    def test_auto_naming(self):
        result = pymafft.align(["ACGT", "ACGT"])
        assert result.sequences[0].name == "seq_1"
        assert result.sequences[1].name == "seq_2"

    def test_identical_sequences(self):
        result = pymafft.align(["ACDEFGHIK", "ACDEFGHIK"])
        assert result.nseq == 2
        assert result.sequences[0].sequence == result.sequences[1].sequence

    def test_residue_preservation(self):
        originals = ["ACDEFGHIKLMNP", "ACDEFHIKLMNP", "ACDEHIKLMNP"]
        result = pymafft.align(originals)
        for i, seq in enumerate(result.sequences):
            ungapped = seq.ungapped()
            assert ungapped == originals[i], (
                f"seq {i}: ungapped '{ungapped}' != original '{originals[i]}'"
            )

    def test_dna_sequences(self):
        result = pymafft.align([
            "ATGCGATCGATCGATCG",
            "ATGCGATCGATCG",
            "ATGCGATCG",
        ])
        assert result.nseq == 3
        assert result.width >= 17  # at least as long as longest

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="no sequences"):
            pymafft.align([])

    def test_invalid_input_type(self):
        with pytest.raises((ValueError, TypeError)):
            pymafft.align("not a list")


class TestStrategies:
    """Test different alignment strategies."""

    def test_fftns2(self):
        result = pymafft.align(["ACDEFGHIK", "ACDEFHIK"], strategy="fftns2")
        assert result.nseq == 2

    def test_fftnsi(self):
        result = pymafft.align(
            ["ACDEFGHIK", "ACDEFHIK", "ACDHIK"],
            strategy="fftnsi",
            maxiterate=2,
        )
        assert result.nseq == 3

    def test_ginsi(self):
        result = pymafft.align(
            ["ACDEFGHIK", "ACDEFHIK", "ACDHIK"],
            strategy="ginsi",
            maxiterate=2,
        )
        assert result.nseq == 3

    def test_linsi(self):
        result = pymafft.align(
            ["ACDEFGHIK", "ACDEFHIK", "ACDHIK"],
            strategy="linsi",
            maxiterate=2,
        )
        assert result.nseq == 3

    def test_einsi(self):
        result = pymafft.align(
            ["ACDEFGHIK", "ACDEFHIK", "ACDHIK"],
            strategy="einsi",
            maxiterate=2,
        )
        assert result.nseq == 3

    def test_invalid_strategy(self):
        with pytest.raises(ValueError, match="unknown strategy"):
            pymafft.align(["ACGT", "ACGT"], strategy="invalid")


class TestAlignFile:
    """Test file-based alignment."""

    @pytest.mark.skipif(
        not os.path.exists(SAMPLE_PATH),
        reason="mafft-upstream submodule not available",
    )
    def test_align_sample(self):
        result = pymafft.align_file(SAMPLE_PATH)
        assert result.nseq == 36
        assert result.width > 0
        for seq in result.sequences:
            assert len(seq) == result.width

    @pytest.mark.skipif(
        not os.path.exists(SAMPLERNA_PATH),
        reason="mafft-upstream submodule not available",
    )
    def test_align_rna(self):
        result = pymafft.align_file(SAMPLERNA_PATH)
        assert result.nseq > 0

    def test_missing_file(self):
        with pytest.raises(ValueError, match="error reading"):
            pymafft.align_file("/nonexistent/file.fa")


class TestAlignFastaString:
    """Test FASTA string alignment."""

    def test_basic(self):
        fasta = ">s1\nACDEFGHIK\n>s2\nACDEFHIK\n>s3\nACDHIK\n"
        result = pymafft.align_fasta_string(fasta)
        assert result.nseq == 3

    def test_preserves_names(self):
        fasta = ">my_protein\nACDEFGHIK\n>other_protein\nACDEFHIK\n"
        result = pymafft.align_fasta_string(fasta)
        assert result.sequences[0].name == "my_protein"
        assert result.sequences[1].name == "other_protein"


# ---------------------------------------------------------------------------
# AlignmentResult API tests
# ---------------------------------------------------------------------------


class TestAlignmentResult:
    """Test AlignmentResult methods."""

    @pytest.fixture
    def result(self):
        return pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])

    def test_len(self, result):
        assert len(result) == 3

    def test_getitem(self, result):
        seq0 = result[0]
        assert isinstance(seq0, pymafft.AlignedSequence)

    def test_getitem_out_of_range(self, result):
        with pytest.raises(ValueError):
            result[99]

    def test_iteration(self, result):
        names = [seq.name for seq in result]
        assert len(names) == 3

    def test_to_fasta(self, result):
        fasta = result.to_fasta()
        assert fasta.count(">") == 3

    def test_to_tuples(self, result):
        tuples = result.to_tuples()
        assert len(tuples) == 3
        assert all(isinstance(t, tuple) and len(t) == 2 for t in tuples)

    def test_repr(self, result):
        r = repr(result)
        assert "AlignmentResult" in r
        assert "nseq=3" in r

    def test_width_property(self, result):
        assert result.width > 0


class TestAlignedSequence:
    """Test AlignedSequence methods."""

    @pytest.fixture
    def seq(self):
        result = pymafft.align(["ACDEFGHIK", "ACDIK"])
        return result[0]

    def test_ungapped(self, seq):
        ungapped = seq.ungapped()
        assert "-" not in ungapped
        assert ungapped == "ACDEFGHIK"

    def test_len(self, seq):
        assert len(seq) > 0

    def test_str(self, seq):
        s = str(seq)
        assert s.startswith(">")

    def test_repr(self, seq):
        r = repr(seq)
        assert "AlignedSequence" in r


# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------


def test_version():
    assert hasattr(pymafft, "__version__")
    assert isinstance(pymafft.__version__, str)
    assert "." in pymafft.__version__
