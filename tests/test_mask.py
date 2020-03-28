import os
import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from augur import mask

# Test inputs for the commands. Writing these here so tests are self-contained.
@pytest.fixture
def sequences():
    return {
        "SEQ1": SeqRecord(Seq("ATGC-ATGC-ATGC"), id="SEQ1"),
        "SEQ2": SeqRecord(Seq("ATATATATATATATAT"), id="SEQ2")
    }

@pytest.fixture
def fasta_file(tmpdir, sequences):
    fasta_file = str(tmpdir / "test.fasta")
    with open(fasta_file, "w") as fh:
        SeqIO.write(sequences.values(), fh, "fasta")
    return fasta_file

TEST_VCF="""\
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
SEQ	1	.	G	A	.		.
SEQ	2	.	G	A	.		.
SEQ	3	.	C	T	.		.
SEQ	5	.	C	T	.		.
SEQ	8	.	A	C	.		.
"""

@pytest.fixture
def vcf_file(tmpdir):
    vcf_file = str(tmpdir / "test.vcf")
    with open(vcf_file, "w") as fh:
        fh.write(TEST_VCF)
    return vcf_file

TEST_BED="""\
Chrom	ChromStart	ChromEnd	locus tag	Comment	
SEQ1	1	2	IG18_Rv0018c-Rv0019c	
SEQ1	4	4	IG18_Rv0018c-Rv0019c	
SEQ1	6	8	IG18_Rv0018c-Rv0019c	
SEQ1	7	10	IG18_Rv0018c-Rv0019c	
"""

@pytest.fixture
def bed_file(tmpdir):
    bed_file = str(tmpdir / "exclude.bed")
    with open(bed_file, "w") as fh:
        fh.write(TEST_BED)
    return bed_file

@pytest.fixture
def mp_context(monkeypatch):
    #Have found explicit monkeypatch context-ing prevents stupid bugs
    with monkeypatch.context() as mp:
        yield mp

class TestMask:
    def test_get_chrom_name_valid(self, vcf_file):
        """get_chrom_name should return the first CHROM field in a vcf file"""
        assert mask.get_chrom_name(vcf_file) == "SEQ"

    def test_get_chrom_name_invalid(self, tmpdir):
        """get_chrom_name should return nothing if no CHROM field is found in the VCF"""
        vcf_file = str(tmpdir / "incomplete.vcf")
        with open(vcf_file, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n")
        assert mask.get_chrom_name(vcf_file) is None
    
    def test_read_bed_file(self, bed_file):
        """read_bed_file should read and deduplicate the list of sites in a bed file"""
        # Not a whole lot of testing to do with bed files. We're basically testing if pandas
        # can read a CSV and numpy can dedupe it.
        assert mask.read_bed_file(bed_file) == [1,2,4,6,7,8,9,10]

    def test_mask_vcf_bails_on_no_chrom(self, tmpdir, mp_context):
        """mask_vcf should pull a sys.exit() if get_chrom_name returns None"""
        bad_vcf = str(tmpdir / "bad.vcf")
        with open(bad_vcf, "w") as fh:
            fh.write("#")
        with pytest.raises(SystemExit) as err:
            mask.mask_vcf([], bad_vcf, "")
    
    def test_mask_vcf_creates_maskfile(self, vcf_file, mp_context):
        """mask_vcf should create and use a mask file from the given list of sites"""
        mask_file = vcf_file + "_maskTemp"
        def shell_has_maskfile(call, **kwargs):
            assert mask_file in call
        mp_context.setattr(mask, "get_chrom_name", lambda f: "SEQ")
        mp_context.setattr(mask, "run_shell_command", shell_has_maskfile)
        mask.mask_vcf([1,5], vcf_file, vcf_file, cleanup=False)
        assert os.path.isfile(mask_file), "Mask file was not written!"
        with open(mask_file) as fh:
            assert fh.read() == "SEQ	1\nSEQ	5", "Incorrect mask file written!"
    
    def test_mask_vcf_handles_gz(self, vcf_file, mp_context):
        """mask_vcf should recognize when the in or out files are .gz and call out accordingly"""
        def test_shell(call, raise_errors=True):
            assert "--gzvcf" in call
            assert "| gzip -c" in call
        mp_context.setattr(mask, "run_shell_command", test_shell)
        mp_context.setattr(mask, "get_chrom_name", lambda f: "SEQ")
        # Still worth using the fixture here because writing the mask file uses this path
        # Using it for both entries in case god forbid someone starts assuming that path
        # valid earlier than it currently is.
        in_file = vcf_file + ".gz"
        mask.mask_vcf([1,5], in_file, in_file)

    def test_mask_vcf_removes_matching_sites(self, tmpdir, vcf_file):
        """mask_vcf should remove the given sites from the VCF file"""
        out_file = str(tmpdir / "output.vcf")
        mask.mask_vcf([5,6], vcf_file, out_file)
        with open(out_file) as after, open(vcf_file) as before:
            assert len(after.readlines()) == len(before.readlines()) - 1, "Too many lines removed!"
            assert "SEQ\t5" not in after.read(), "Correct sites not removed!"

    def test_mask_vcf_cleanup_flag(self, vcf_file, mp_context):
        """mask_vcf should respect the cleanup flag"""
        tmp_mask_file = vcf_file + "_maskTemp"
        mp_context.setattr(mask, "run_shell_command", lambda *a, **k: None)
        mp_context.setattr(mask, "get_chrom_name", lambda f: "SEQ")

        mask.mask_vcf([], vcf_file, "", cleanup=True)
        assert not os.path.isfile(tmp_mask_file), "Temporary mask not cleaned up"

        mask.mask_vcf([], vcf_file, "", cleanup=False)
        assert os.path.isfile(tmp_mask_file), "Temporary mask cleaned up as expected"
    
    def test_mask_fasta_normal_case(self, tmpdir, fasta_file, sequences):
        """mask_fasta normal case - all sites in sequences"""
        out_file = str(tmpdir / "output.fasta")
        mask_sites = [5,10]
        mask.mask_fasta([5,10], fasta_file, out_file)
        output = SeqIO.parse(out_file, "fasta")
        for seq in output:
            original = sequences[seq.id]
            for idx, site in enumerate(seq):
                if idx not in mask_sites:
                    assert site == original[idx], "Incorrect sites modified!"
                else:
                    assert site == "N", "Not all sites modified correctly!"
    
    def test_mask_fasta_out_of_index(self, tmpdir, fasta_file, sequences):
        """mask_fasta provided a list of indexes past the length of the sequences"""
        out_file = str(tmpdir / "output.fasta")
        max_length = max(len(record.seq) for record in sequences.values())
        mask.mask_fasta([5, max_length, max_length+5], fasta_file, out_file)
        output = SeqIO.parse(out_file, "fasta")
        for seq in output:
            assert seq[5] == "N", "Not all sites masked correctly!"
            original = sequences[seq.id]
            for idx, site in enumerate(seq):
                if idx != 5:
                    assert site == original[idx], "Incorrect sites modified!"