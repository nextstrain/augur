from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from shlex import quote
import pathlib

from augur import align

import pytest


class TestAlign:
    def test_make_gaps_ambiguous(self):
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("G-AC")), SeqRecord(Seq("----")), SeqRecord(Seq("TAGC"))]
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
            
    def test_prune_seqs_matching_alignment(self):
        sequence = {
            "seq1": SeqRecord(Seq("GTAC"), name="seq1"),
            "seq2": SeqRecord(Seq("CGTT"), name="seq2"),
            "seq3": SeqRecord(Seq("TAGC"), name="seq3"),
        }
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("GTAC"), name="seq1"),
                SeqRecord(Seq("TAGC"), name="seq3"),
            ]
        )
        
        result = align.prune_seqs_matching_alignment(sequence, alignment)
        assert list(result.keys()) == ["seq2"]
        assert result["seq2"].seq == sequence["seq2"].seq

    def test_generate_alignment_cmd_non_mafft(self):
        with pytest.raises(align.AlignmentError):
            assert align.generate_alignment_cmd('no-mafft', 1, None, None, None, None)
            
    def test_generate_alignment_cmd_mafft_existing_aln_fname(self):
        existing_aln_fname = "existing_aln"
        seqs_to_align_fname = "seqs_to_align"
        aln_fname = "aln_fname"
        log_fname = "log_fname"
        
        result = align.generate_alignment_cmd("mafft", 1,
                                              existing_aln_fname,
                                              seqs_to_align_fname,
                                              aln_fname,
                                              log_fname)
        
        expected = "mafft --add %s --keeplength --reorder --anysymbol --nomemsave --thread %d %s 1> %s 2> %s" % (quote(seqs_to_align_fname), 1, quote(existing_aln_fname), quote(aln_fname), quote(log_fname))
        
        assert result == expected
                                    
    def test_generate_alignment_cmd_mafft_no_existing_aln_fname(self):
        seqs_to_align_fname = "seqs_to_align"
        aln_fname = "aln_fname"
        log_fname = "log_fname"
        
        result = align.generate_alignment_cmd("mafft", 1,
                                              None,
                                              seqs_to_align_fname,
                                              aln_fname,
                                              log_fname)
        
        expected = "mafft --reorder --anysymbol --nomemsave --thread %d %s 1> %s 2> %s" % (1, quote(seqs_to_align_fname), quote(aln_fname), quote(log_fname))
        
        assert result == expected
        
    def test_read_alignment(self):
        data_file = pathlib.Path('tests/data/align/test_aligned_sequences.fasta')
        result = align.read_alignment(str(data_file.resolve()))
        
        assert len(result) == 3
        
    def test_read_sequences(self):
        data_file = pathlib.Path('tests/data/align/test_aligned_sequences.fasta')
        result = align.read_sequences(data_file)
        assert len(result.keys()) == 3
        
        
