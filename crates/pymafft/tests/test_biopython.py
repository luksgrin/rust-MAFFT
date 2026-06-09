"""Biopython interop tests.

pymafft duck-types `Bio.SeqRecord.SeqRecord` (and anything else with `.id` /
`.seq` attributes), so the typical Biopython workflow works without extra
glue:

    records = list(SeqIO.parse("input.fasta", "fasta"))
    result = pymafft.align(records, strategy="linsi", maxiterate=1000)

The output side is FASTA-string-based: convert via `SeqIO.parse(StringIO(
result.to_fasta()), "fasta")` to get a list of `SeqRecord` back, or build a
`MultipleSeqAlignment` directly from the result.

The whole suite skips if biopython isn't installed (it's an optional dep).
"""

from __future__ import annotations

from io import StringIO

import pytest


biopython = pytest.importorskip("Bio")
import pymafft  # noqa: E402  (after the skip)
from Bio import AlignIO, SeqIO  # noqa: E402
from Bio.Align import MultipleSeqAlignment  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def records():
    """Three protein SeqRecords (alpha, beta, gamma)."""
    return [
        SeqRecord(Seq("ACDEFGHIK"), id="alpha", description="protein alpha"),
        SeqRecord(Seq("ACDEFHIK"), id="beta", description="protein beta"),
        SeqRecord(Seq("ACDHIK"), id="gamma", description="protein gamma"),
    ]


# ---------------------------------------------------------------------------
# Direct SeqRecord input
# ---------------------------------------------------------------------------


class TestSeqRecordInput:
    """pymafft accepts a list of `SeqRecord` directly."""

    def test_align_accepts_seqrecord_list(self, records):
        result = pymafft.align(records)
        assert result.nseq == 3
        assert [s.name for s in result.sequences] == ["alpha", "beta", "gamma"]

    def test_align_with_strategy(self, records):
        result = pymafft.align(records, strategy="linsi", maxiterate=1000)
        assert result.nseq == 3
        # ungapped residues equal original .seq
        for r, s in zip(records, result.sequences):
            assert s.ungapped() == str(r.seq)

    def test_align_with_all_strategies(self, records):
        for strategy in ["fftns2", "fftnsi", "ginsi", "linsi", "einsi"]:
            result = pymafft.align(records, strategy=strategy)
            assert result.nseq == 3, f"{strategy} dropped sequences"

    def test_seqrecord_via_iter(self, records):
        """Iterators, not just lists, are accepted."""
        result = pymafft.align(iter(records))
        assert result.nseq == 3

    def test_mixed_input_rejected(self, records):
        """Mixing SeqRecord with bare strings is rejected (consistent shape)."""
        bad = [records[0], "ACGT"]
        with pytest.raises(ValueError, match=r"is not a string"):
            pymafft.align(bad)

    def test_objects_without_seq_attr_rejected(self):
        """Plain objects without `.id`/`.seq` raise a clear error."""
        class NotASeqRecord:
            pass

        with pytest.raises(ValueError, match=r"is not a string"):
            pymafft.align([NotASeqRecord(), NotASeqRecord()])

    def test_id_taken_from_record(self, records):
        """`SeqRecord.id` becomes the alignment name (not `.name` or `.description`)."""
        records[0].id = "renamed_alpha"
        records[0].name = "different_name"
        result = pymafft.align(records)
        assert result.sequences[0].name == "renamed_alpha"


# ---------------------------------------------------------------------------
# Output → Biopython
# ---------------------------------------------------------------------------


class TestOutputToBiopython:
    """Converting `AlignmentResult` into Biopython types."""

    def test_to_fasta_parses_with_seqio(self, records):
        result = pymafft.align(records)
        fasta_text = result.to_fasta()
        parsed = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
        assert len(parsed) == 3
        assert all(isinstance(r, SeqRecord) for r in parsed)
        assert [r.id for r in parsed] == ["alpha", "beta", "gamma"]

    def test_to_fasta_parses_as_alignio(self, records):
        """`AlignIO.read` returns a `MultipleSeqAlignment`."""
        result = pymafft.align(records)
        msa = AlignIO.read(StringIO(result.to_fasta()), "fasta")
        assert isinstance(msa, MultipleSeqAlignment)
        assert len(msa) == 3
        assert msa.get_alignment_length() == result.width

    def test_manual_msa_construction(self, records):
        """Build a `MultipleSeqAlignment` directly from `.sequences`."""
        result = pymafft.align(records)
        aligned_records = [
            SeqRecord(Seq(s.sequence), id=s.name, description="")
            for s in result.sequences
        ]
        msa = MultipleSeqAlignment(aligned_records)
        assert len(msa) == 3
        assert msa.get_alignment_length() == result.width

    def test_to_tuples_dict_round_trip(self, records):
        """`dict(result.to_tuples())` gives a {name: aligned_seq} map."""
        result = pymafft.align(records)
        m = dict(result.to_tuples())
        assert set(m) == {"alpha", "beta", "gamma"}
        # The unaligned-residue content still matches the original input.
        for r in records:
            assert m[r.id].replace("-", "") == str(r.seq)


# ---------------------------------------------------------------------------
# Round-trip parity
# ---------------------------------------------------------------------------


class TestBiopythonRoundTrip:
    """Equivalent results across input-conversion paths."""

    def test_seqrecord_vs_tuple_paths_agree(self, records):
        """SeqRecord input and tuple input give the same alignment."""
        from_records = pymafft.align(records)
        from_tuples = pymafft.align([(r.id, str(r.seq)) for r in records])
        assert from_records.to_tuples() == from_tuples.to_tuples()

    def test_seqrecord_vs_fasta_string_paths_agree(self, records):
        """SeqRecord input and FASTA-string input give the same aligned residues.

        Note: `SeqIO.write` emits FASTA headers as `id<space>description`,
        so the resulting NAMES differ from the SeqRecord-input path (which
        uses just `id`). Compare aligned sequences only.
        """
        from_records = pymafft.align(records)
        buf = StringIO()
        SeqIO.write(records, buf, "fasta")
        from_string = pymafft.align_fasta_string(buf.getvalue())
        assert (
            [s for _, s in from_records.to_tuples()]
            == [s for _, s in from_string.to_tuples()]
        )

    def test_fasta_file_via_seqio_matches_align_file(self, tmp_path, records):
        """`SeqIO.write` + `align_file` agrees with `align(records)` on residues.

        Same caveat as `test_seqrecord_vs_fasta_string_paths_agree`: names
        come from different sources, so compare sequences only.
        """
        path = tmp_path / "input.fasta"
        SeqIO.write(records, str(path), "fasta")
        from_file = pymafft.align_file(str(path))
        from_records = pymafft.align(records)
        assert (
            [s for _, s in from_file.to_tuples()]
            == [s for _, s in from_records.to_tuples()]
        )

    def test_msa_to_records_to_align(self, records):
        """`MultipleSeqAlignment` → records → `pymafft.align`."""
        # An MSA needs same-length records; pad to alignment width.
        result = pymafft.align(records)
        aligned_records = [
            SeqRecord(Seq(s.sequence), id=s.name, description="")
            for s in result.sequences
        ]
        msa = MultipleSeqAlignment(aligned_records)
        # Feeding the MSA back into pymafft is the same as re-aligning the
        # already-aligned records (and dropping the gaps first).
        unaligned = [
            SeqRecord(Seq(str(r.seq).replace("-", "")), id=r.id)
            for r in msa
        ]
        result2 = pymafft.align(unaligned)
        # Ungapped residues identical to original.
        for original, after in zip(records, result2.sequences):
            assert after.ungapped() == str(original.seq)


# ---------------------------------------------------------------------------
# Workflow examples (smoke tests for documented patterns)
# ---------------------------------------------------------------------------


class TestDocumentedWorkflows:
    """Smoke tests for the workflows the README will document."""

    def test_workflow_read_align_write(self, tmp_path):
        """SeqIO read → pymafft align → SeqIO write."""
        # Set up: write a small input file
        records_in = [
            SeqRecord(Seq("MNGTEGDNFYVPFSNKTG"), id="opsin1"),
            SeqRecord(Seq("MNGTEGPNFYVPFSNITG"), id="opsin2"),
            SeqRecord(Seq("MNGTEGINFYVPMSNKTG"), id="opsin3"),
        ]
        input_path = tmp_path / "opsins.fasta"
        SeqIO.write(records_in, str(input_path), "fasta")

        # Workflow
        records = list(SeqIO.parse(str(input_path), "fasta"))
        result = pymafft.align(records, strategy="linsi", maxiterate=1000)
        aligned_records = [
            SeqRecord(Seq(s.sequence), id=s.name)
            for s in result.sequences
        ]
        output_path = tmp_path / "opsins.aligned.fasta"
        SeqIO.write(aligned_records, str(output_path), "fasta")

        # Verify the round-trip
        reread = list(SeqIO.parse(str(output_path), "fasta"))
        assert len(reread) == 3
        widths = {len(r.seq) for r in reread}
        assert len(widths) == 1, "all aligned rows must have same length"

    def test_workflow_direct_to_alignio(self, records):
        """Alternative: skip the in-memory record list, go via FASTA string."""
        result = pymafft.align(records)
        msa = AlignIO.read(StringIO(result.to_fasta()), "fasta")
        assert isinstance(msa, MultipleSeqAlignment)
        # All Biopython column ops work on the MSA.
        assert msa[:, 0] == "AAA"  # first column: all three records start with 'A'
