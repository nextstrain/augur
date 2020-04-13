import os

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from shlex import quote

from augur import align

import pytest
import pathlib

def write_strains(tmpdir, name, strains):
    path = str(tmpdir / name + ".fasta")
    with open(path, "w") as fh:
        SeqIO.write(strains, fh, "fasta")
    return path

@pytest.fixture
def ref_seq():
    return SeqRecord(Seq("AAAATTTTGGGGCCCC"), "REF")

@pytest.fixture
def test_seqs(ref_seq):
    return {
        "PREFIX": SeqRecord(ref_seq.seq[3:], "PREFIX"),
        "SUFFIX": SeqRecord(ref_seq.seq[:-3], "SUFFIX"),
        "LONGER": SeqRecord("CCC" + ref_seq.seq + "AAA", "LONGER")
    }

@pytest.fixture
def existing_aln(ref_seq):
    return {
        "EXISTING1": SeqRecord("NN" + ref_seq.seq[2:-3] + "NNN", "EXISTING1"),
        "EXISTING2": SeqRecord("NNN" + ref_seq.seq[3:-1] + "N", "EXISTING2")
    }

@pytest.fixture
def ref_file(tmpdir, ref_seq):
    return write_strains(tmpdir, "ref", ref_seq)

@pytest.fixture
def test_file(tmpdir, test_seqs):
    return write_strains(tmpdir, "test", test_seqs.values())

@pytest.fixture
def test_with_ref(tmpdir, test_seqs, ref_seq):
    return write_strains(tmpdir, "test_w_ref", [ref_seq,] + list(existing_aln.values()))

@pytest.fixture
def existing_file(tmpdir, existing_aln):
    return write_strains(tmpdir, "existing", existing_aln.values())

@pytest.fixture
def existing_with_ref(tmpdir, existing_aln, ref_seq):
    return write_strains(tmpdir, "existing_w_ref", [ref_seq,] + list(existing_aln.values()))

@pytest.fixture
def out_file(tmpdir):
    out_file = str(tmpdir / "out")
    open(out_file, "w").close()
    return out_file

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
        
        result = align.prune_seqs_matching_alignment(sequence.values(), alignment)
        assert [r.name for r in result] == ["seq2"]
        for r in result:
            assert r.seq == sequence[r.name].seq

    def test_prettify_alignment(self):
        data_file = pathlib.Path('tests/data/align/test_aligned_sequences.fasta')
        alignment = align.read_alignment(str(data_file.resolve()))
        seqs = {s.id:s for s in alignment}
        assert "_R_crick_strand" in seqs

        align.prettify_alignment(alignment)
        seqs = {s.id:s for s in alignment}
        assert "crick_strand" in seqs

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
        
        expected = "mafft --add %s --keeplength --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s" % (quote(seqs_to_align_fname), 1, quote(existing_aln_fname), quote(aln_fname), quote(log_fname))
        
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
        
        expected = "mafft --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s" % (1, quote(seqs_to_align_fname), quote(aln_fname), quote(log_fname))
        
        assert result == expected
        
    def test_read_alignment(self):
        data_file = pathlib.Path('tests/data/align/test_aligned_sequences.fasta')
        result = align.read_alignment(str(data_file.resolve()))
        
        assert len(result) == 4
        
    def test_read_sequences(self):
        data_file = pathlib.Path('tests/data/align/test_aligned_sequences.fasta')
        result = align.read_sequences(data_file)
        assert len(result) == 4

    def test_read_seq_compare(self):
        data_file = pathlib.Path("tests/data/align/aa-seq_h3n2_ha_2y_2HA1_dup.fasta")
        with pytest.raises(align.AlignmentError):
            assert align.read_sequences(data_file)

    def test_prepare_no_alignment_or_ref(self, test_file, test_seqs, tmpdir, out_file):
        """
        Strictly we shouldn't see this case, since we should always have a
        ref_name or ref_seq, but it's still good to check the behavior.
        """
        expected_output = out_file + ".to_align.fasta"
        align.prepare([test_file,], None, out_file, None, None)
        assert os.path.isfile(expected_output), "Didn't write sequences where we expected"
        for name, seq in SeqIO.to_dict(SeqIO.parse(expected_output, "fasta")).items():
            assert seq.seq == test_seqs[name].seq