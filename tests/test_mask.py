"""Tests for augur mask

NOTE: Several functions are monkeypatched in these tests. If you change the arguments
for any function in mask.py, check that it is correctly updated in this file.
"""
import argparse
import os
import pytest
import string

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from augur import mask
from augur.utils import VALID_NUCLEOTIDES

# Test inputs for the commands. Writing these here so tests are self-contained.
@pytest.fixture
def sequences():
    return {
        "SEQ1": SeqRecord(Seq("ATGC-ATGC-ATGC"), id="SEQ1"),
        "SEQ2": SeqRecord(Seq("ATATATATATATATAT"), id="SEQ2"),
        "SEQ_ALLCHARS": SeqRecord(Seq(string.ascii_letters + string.digits + "-"), id="SEQ_ALLCHARS")
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
SEQ	13	.	A	C	.		.
"""

@pytest.fixture
def vcf_file(tmpdir):
    vcf_file = str(tmpdir / "test.vcf")
    with open(vcf_file, "w") as fh:
        fh.write(TEST_VCF)
    return vcf_file

TEST_BED_SEQUENCE = [1,4,6,7,8,9]
# IF YOU UPDATE ONE OF THESE, UPDATE THE OTHER.
TEST_BED="""\
SEQ1	1	2	IG18_Rv0018c-Rv0019c	
SEQ1	4	5	IG18_Rv0018c-Rv0019c	
SEQ1	6	8	IG18_Rv0018c-Rv0019c	
SEQ1	7	10	IG18_Rv0018c-Rv0019c	
"""

@pytest.fixture
def bed_file(tmpdir):
    bed_file = str(tmpdir / "exclude.bed")
    with open(bed_file, "w") as fh:
        fh.write(TEST_BED)
    return bed_file

TEST_MASK_SEQUENCE = [3,5,7]
TEST_MASK="""\
4
6
8
"""

@pytest.fixture
def mask_file(tmpdir):
    mask_file = str(tmpdir / "exclude.mask")
    with open(mask_file, "w") as fh:
        fh.write(TEST_MASK)
    return mask_file

@pytest.fixture
def out_file(tmpdir):
    out_file = str(tmpdir / "out")
    open(out_file, "w").close()
    return out_file

@pytest.fixture
def mp_context(monkeypatch):
    #Have found explicit monkeypatch context-ing prevents stupid bugs
    with monkeypatch.context() as mp:
        yield mp


@pytest.fixture
def argparser():
    """Provide an easy way to test command line arguments"""
    parser = argparse.ArgumentParser()
    mask.register_arguments(parser)
    def parse(args):
        return parser.parse_args(args.split(" "))
    return parse

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
    
    def test_mask_vcf_bails_on_no_chrom(self, tmpdir):
        """mask_vcf should pull a sys.exit() if get_chrom_name returns None"""
        bad_vcf = str(tmpdir / "bad.vcf")
        with open(bad_vcf, "w") as fh:
            fh.write("#")
        with pytest.raises(SystemExit) as err:
            mask.mask_vcf([], bad_vcf, "")
    
    def test_mask_vcf_creates_maskfile(self, vcf_file, mp_context):
        """mask_vcf should create a 1-indexed mask file from the 0-indexed list of sites"""
        mask_file = vcf_file + "_maskTemp"
        def shell_has_maskfile(call, **kwargs):
            assert mask_file in call
        mp_context.setattr(mask, "get_chrom_name", lambda f: "SEQ")
        mp_context.setattr(mask, "run_shell_command", shell_has_maskfile)
        mask.mask_vcf([1,5], vcf_file, vcf_file, cleanup=False)
        assert os.path.isfile(mask_file), "Mask file was not written!"
        with open(mask_file) as fh:
            assert fh.read() == "SEQ	2\nSEQ	6", "Incorrect mask file written!"
    
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

    def test_mask_vcf_removes_matching_sites(self, vcf_file, out_file):
        """mask_vcf should remove the given sites from the VCF file"""
        mask.mask_vcf([4,5], vcf_file, out_file)
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
    
    def test_mask_fasta_normal_case(self, fasta_file, out_file, sequences):
        """mask_fasta normal case - all sites in sequences"""
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
    
    def test_mask_fasta_out_of_index(self, out_file, fasta_file, sequences):
        """mask_fasta provided a list of indexes past the length of the sequences"""
        max_length = max(len(record.seq) for record in sequences.values())
        mask.mask_fasta([5, max_length, max_length+5], fasta_file, out_file)
        output = SeqIO.parse(out_file, "fasta")
        for seq in output:
            assert seq[5] == "N", "Not all sites masked correctly!"
            original = sequences[seq.id]
            for idx, site in enumerate(seq):
                if idx != 5:
                    assert site == original[idx], "Incorrect sites modified!"
    
    def test_mask_fasta_from_beginning(self, out_file, fasta_file, sequences):
        mask.mask_fasta([], fasta_file, out_file, mask_from_beginning=3)
        output = SeqIO.parse(out_file, "fasta")
        for seq in output:
            original = sequences[seq.id]
            assert seq.seq[:3] == "NNN"
            assert seq.seq[3:] == original.seq[3:]

    def test_mask_fasta_from_end(self, out_file, fasta_file, sequences):
        mask.mask_fasta([], fasta_file, out_file, mask_from_end=3)
        output = SeqIO.parse(out_file, "fasta")
        for seq in output:
            original = sequences[seq.id]
            assert seq.seq[-3:] == "NNN"
            assert seq.seq[:-3] == original.seq[:-3]

    def test_mask_fasta_from_beginning_and_end(self, out_file, fasta_file, sequences):
        mask.mask_fasta([], fasta_file, out_file, mask_from_beginning=2, mask_from_end=3)
        output = SeqIO.parse(out_file, "fasta")
        for seq in output:
            original = sequences[seq.id]
            assert seq.seq[:2] == "NN"
            assert seq.seq[-3:] == "NNN"
            assert seq.seq[2:-3] == original.seq[2:-3]

    @pytest.mark.parametrize("beginning,end", ((1000,0), (0,1000),(1000,1000)))
    def test_mask_fasta_from_beginning_and_end_too_long(self, fasta_file, out_file, beginning, end):
        mask.mask_fasta([], fasta_file, out_file, mask_from_beginning=beginning, mask_from_end=end)
        output = SeqIO.parse(out_file, "fasta")
        for record in output:
            assert record.seq == "N" * len(record.seq)

    @pytest.mark.parametrize("mask_invalid", (True, False))
    def test_mask_fasta_invalid_sites(self, fasta_file, out_file, sequences, mask_invalid):
        """Verify that mask_fasta masks invalid nucleotides when and only when mask_invalid is passed as True"""
        mask.mask_fasta([], fasta_file, out_file, mask_invalid=mask_invalid)
        output = SeqIO.to_dict(SeqIO.parse(out_file, "fasta"))
        for name, record in sequences.items():
            for site, nucleotide in enumerate(record.seq):
                if nucleotide not in VALID_NUCLEOTIDES and mask_invalid is True:
                    assert output[name][site] == "N"
                else:
                    assert output[name][site] == nucleotide

    def test_run_handle_missing_sequence_file(self, vcf_file, argparser):
        os.remove(vcf_file)
        args = argparser("-s %s" % vcf_file)
        with pytest.raises(SystemExit):
            mask.run(args)

    def test_run_handle_empty_sequence_file(self, vcf_file, argparser):
        open(vcf_file,"w").close()
        args = argparser("-s %s --mask-sites 1" % vcf_file)
        with pytest.raises(SystemExit):
            mask.run(args)

    def test_run_handle_missing_mask_file(self, vcf_file, bed_file, argparser):
        os.remove(bed_file)
        args = argparser("-s %s --mask %s" % (vcf_file, bed_file))
        with pytest.raises(SystemExit):
            mask.run(args)

    def test_run_handle_empty_mask_file(self, vcf_file, bed_file, argparser):
        open(bed_file, "w").close()
        args = argparser("-s %s --mask %s" % (vcf_file, bed_file))
        with pytest.raises(SystemExit):
            mask.run(args)

    def test_run_recognize_vcf(self, bed_file, vcf_file, argparser, mp_context):
        """Ensure we're handling vcf files correctly"""
        args = argparser("--mask=%s -s %s --no-cleanup" % (bed_file, vcf_file))
        def fail(*args, **kwargs):
            assert False, "Called mask_fasta incorrectly"
        mp_context.setattr(mask, "mask_vcf", lambda *a, **k: None)
        mp_context.setattr(mask, "mask_fasta", fail)
        mp_context.setattr(mask, "copyfile", lambda *args: None)
        mask.run(args)

    def test_run_recognize_fasta(self, bed_file, fasta_file, argparser, mp_context):
        """Ensure we're handling fasta files correctly"""
        args = argparser("--mask=%s -s %s --no-cleanup" % (bed_file, fasta_file))
        def fail(*args, **kwargs):
            assert False, "Called mask_fasta incorrectly"
        mp_context.setattr(mask, "mask_fasta", lambda *a, **k: None)
        mp_context.setattr(mask, "mask_vcf", fail)
        mp_context.setattr(mask, "copyfile", lambda *args: None)
        mask.run(args)

    def test_run_handle_missing_outfile(self, bed_file, fasta_file, argparser, mp_context):
        args = argparser("--mask=%s -s %s" % (bed_file, fasta_file))
        expected_outfile = os.path.join(os.path.dirname(fasta_file), "masked_" + os.path.basename(fasta_file))
        def check_outfile(mask_sites, in_file, out_file, **kwargs):
            assert out_file == expected_outfile
            with open(out_file, "w") as fh:
                fh.write("test_string")
        mp_context.setattr(mask, "mask_fasta", check_outfile)
        mask.run(args)
        with open(fasta_file) as fh:
            assert fh.read() == "test_string"
    
    def test_run_respect_no_cleanup(self, bed_file, vcf_file, argparser, mp_context):
        out_file = os.path.join(os.path.dirname(vcf_file), "masked_" + os.path.basename(vcf_file))
        def make_outfile(mask_sites, in_file, out_file, cleanup=True):
            assert cleanup == False
            open(out_file, "w").close() # need out_file to exist
        mp_context.setattr(mask, "mask_vcf", make_outfile)
        args = argparser("--mask=%s -s %s -o %s --no-cleanup" % (bed_file, vcf_file, out_file))
        mask.run(args)
        assert os.path.exists(out_file), "Output file incorrectly deleted"

    def test_run_normal_case(self, bed_file, vcf_file, out_file, argparser, mp_context):
        def check_args(mask_sites, in_file, _out_file, cleanup):
            assert mask_sites == TEST_BED_SEQUENCE, "Wrong mask sites provided"
            assert in_file == vcf_file, "Incorrect input file provided"
            assert _out_file == out_file, "Incorrect output file provided"
            assert cleanup is True, "Cleanup erroneously passed in as False"
        mp_context.setattr(mask, "mask_vcf", check_args)
        args = argparser("--mask=%s --sequences=%s --output=%s" %(bed_file, vcf_file, out_file))
        mask.run(args)
        assert os.path.exists(out_file), "Output file incorrectly deleted"
    
    def test_run_with_mask_sites(self, vcf_file, out_file, argparser, mp_context):
        args = argparser("--mask-sites 2 8 -s %s -o %s" % (vcf_file, out_file))
        def check_mask_sites(mask_sites, *args, **kwargs):
            # mask-sites are passed to the CLI as one-indexed
            assert mask_sites == [1,7]
        mp_context.setattr(mask, "mask_vcf", check_mask_sites)
        mask.run(args)

    def test_run_with_mask_sites_and_mask_file(self, vcf_file, out_file, bed_file, argparser, mp_context):
        args = argparser("--mask-sites 20 21 --mask %s -s %s -o %s" % (bed_file, vcf_file, out_file))
        def check_mask_sites(mask_sites, *args, **kwargs):
            # mask-sites are passed to the CLI as one-indexed
            assert mask_sites == sorted(set(TEST_BED_SEQUENCE + [19,20]))
        mp_context.setattr(mask, "mask_vcf", check_mask_sites)
        mask.run(args)

    def test_run_requires_some_masking(self, vcf_file, argparser):
        args = argparser("-s %s" % vcf_file)
        with pytest.raises(SystemExit) as err:
            mask.run(args)

    @pytest.mark.parametrize("op", ("beginning", "end"))
    def test_run_vcf_cannot_mask_beginning_or_end(self, vcf_file, argparser, op):
        args = argparser("-s %s --mask-from-%s 2" % (vcf_file, op))
        with pytest.raises(SystemExit) as err:
            mask.run(args)

    def test_run_fasta_mask_from_beginning_or_end(self, fasta_file, out_file, argparser, mp_context):
        args = argparser("-s %s -o %s --mask-from-beginning 2 --mask-from-end 3" % (fasta_file, out_file))
        def check_mask_from(*args, mask_from_beginning=0, mask_from_end=0, **kwargs):
            assert mask_from_beginning == 2
            assert mask_from_end == 3
        mp_context.setattr(mask, "mask_fasta", check_mask_from)
        mask.run(args)

    def test_e2e_fasta_minimal(self, fasta_file, bed_file, sequences, argparser):
        args = argparser("-s %s --mask %s" % (fasta_file, bed_file))
        mask.run(args)
        output = SeqIO.parse(fasta_file,"fasta")
        for record in output:
            reference = sequences[record.id].seq
            for idx, site in enumerate(record.seq):
                if idx in TEST_BED_SEQUENCE:
                    assert site == "N"
                else:
                    assert site == reference[idx]

    def test_e2e_fasta_mask_file(self, fasta_file, mask_file, sequences, argparser):
        args = argparser("-s %s --mask %s" % (fasta_file, mask_file))
        mask.run(args)
        output = SeqIO.parse(fasta_file,"fasta")
        for record in output:
            reference = sequences[record.id].seq
            for idx, site in enumerate(record.seq):
                if idx in TEST_MASK_SEQUENCE:
                    assert site == "N"
                else:
                    assert site == reference[idx]

    def test_e2e_fasta_beginning_end_sites(self, fasta_file, bed_file, out_file, sequences, argparser):
        from_beginning = 3
        from_end = 1
        arg_sites = [5, 12]
        expected_removals = sorted(set(TEST_BED_SEQUENCE + [s - 1 for s in arg_sites]))
        print(expected_removals)
        args = argparser("-s %s -o %s --mask %s --mask-from-beginning %s --mask-from-end %s --mask-sites %s" % (
                         fasta_file, out_file, bed_file, from_beginning, from_end, " ".join(str(s) for s in arg_sites)))
        mask.run(args)
        output = SeqIO.parse(out_file, "fasta")
        for record in output:
            reference = str(sequences[record.id].seq)
            masked_seq = str(record.seq)
            assert masked_seq[:from_beginning] == "N" * from_beginning
            assert masked_seq[-from_end:] == "N" * from_end
            for idx, site in enumerate(masked_seq[from_beginning:-from_end], from_beginning):
                if idx in expected_removals:
                    assert site == "N"
                else:
                    assert site == reference[idx]

    def test_e2e_fasta_mask_invalid(self, fasta_file, out_file, sequences, argparser):
        args = argparser("-s %s -o %s --mask-invalid" % (fasta_file, out_file))
        mask.run(args)
        output = SeqIO.parse(out_file, "fasta")
        for record in output:
            reference = str(sequences[record.id].seq)
            for idx, site in enumerate(reference):
                assert record.seq[idx] == site if site in VALID_NUCLEOTIDES else "N"

    def test_e2e_vcf_minimal(self, vcf_file, bed_file, argparser):
        args = argparser("-s %s --mask %s" % (vcf_file, bed_file))
        mask.run(args)
        with open(vcf_file) as output:
            assert output.readline().startswith("##fileformat") # is a VCF
            assert output.readline().startswith("#CHROM\tPOS\t") # have a header
            for line in output.readlines():
                site = int(line.split("\t")[1]) # POS column
                site = site - 1 # shift to zero-indexed site
                assert site not in TEST_BED_SEQUENCE

    def test_e2e_vcf_mask_file(self, vcf_file, mask_file, argparser):
        args = argparser("-s %s --mask %s" % (vcf_file, mask_file))
        mask.run(args)
        with open(vcf_file) as output:
            assert output.readline().startswith("##fileformat") # is a VCF
            assert output.readline().startswith("#CHROM\tPOS\t") # have a header
            for line in output.readlines():
                site = int(line.split("\t")[1]) # POS column
                site = site - 1 # shift to zero-indexed site
                assert site not in TEST_MASK_SEQUENCE

    def test_e2e_vcf_with_options(self, vcf_file, bed_file, out_file, argparser):
        arg_sites = [5, 12, 14]
        expected_removals = sorted(set(TEST_BED_SEQUENCE + [s - 1 for s in arg_sites]))
        args = argparser("-s %s -o %s --mask %s --mask-sites %s" % (
                         vcf_file, out_file, bed_file, " ".join(str(s) for s in arg_sites)))
        mask.run(args)
        with open(out_file) as output:
            assert output.readline().startswith("##fileformat") # is a VCF
            assert output.readline().startswith("#CHROM\tPOS\t") # have a header
            for line in output.readlines():
                site = int(line.split("\t")[1]) # POS column
                site = site - 1 #re-zero-index the VCF sites
                assert site not in expected_removals
