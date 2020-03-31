from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from augur import align

import pytest


class TestAlign:
    def test_make_gaps_ambiguous(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("G-AC")), SeqRecord(Seq("----")),
                SeqRecord(Seq("TAGC"))
            ]
        )

        align.make_gaps_ambiguous(alignment)

        assert [seq_record.seq for seq_record in alignment] == [
            Seq("GNAC"),
            Seq("NNNN"),
            Seq("TAGC"),
        ]

    def test_check_duplicates_no_arguments(self):
        assert align.check_duplicates() is None

    def test_check_duplicates_strings_with_no_duplicates(self):
        assert align.check_duplicates("GTAC", "CGTT") is None

    def test_check_duplicates_MSA_with_no_duplicates(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("GTAC"), name="seq1"),
                SeqRecord(Seq("CGTT"), name="seq2"),
                SeqRecord(Seq("TAGC"), name="seq3"),
            ]
        )
        assert align.check_duplicates(alignment) is None

    def test_check_duplicates_MSA_and_string_with_no_duplicates(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("GTAC"), name="seq1"),
                SeqRecord(Seq("CGTT"), name="seq2"),
                SeqRecord(Seq("TAGC"), name="seq3"),
            ]
        )
        assert align.check_duplicates(alignment, "TGTT") is None

    def test_check_duplicates_string_with_duplicates(self):
        with pytest.raises(align.AlignmentError):
            assert align.check_duplicates("GTAC", "CGTT", "CGTT")

    def test_check_duplicates_MSA_with_duplicates(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("GTAC"), name="seq1"),
                SeqRecord(Seq("CGTT"), name="seq2"),
                SeqRecord(Seq("TAGC"), name="seq3"),
                SeqRecord(Seq("TAGC"), name="seq3"),
            ]
        )
        with pytest.raises(align.AlignmentError):
            assert align.check_duplicates(alignment)

    def test_check_duplicates_MSA_and_string_with_duplicates(self):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("GTAC"), name="seq1"),
                SeqRecord(Seq("CGTT"), name="seq2"),
                SeqRecord(Seq("TAGC"), name="seq3"),
            ]
        )
        with pytest.raises(align.AlignmentError):
            assert align.check_duplicates(alignment, "seq3")
