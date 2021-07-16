import argparse
import random
import shlex

import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import augur.filter
from augur.utils import read_metadata

@pytest.fixture
def argparser():
    parser = argparse.ArgumentParser()
    augur.filter.register_arguments(parser)
    def parse(args):
        return parser.parse_args(shlex.split(args))
    return parse

@pytest.fixture
def sequences():
    def random_seq(k):
        return "".join(random.choices(("A","T","G","C"), k=k))
    return {
        "SEQ_1": SeqRecord(Seq(random_seq(10)), id="SEQ_1"),
        "SEQ_2": SeqRecord(Seq(random_seq(10)), id="SEQ_2"),
        "SEQ_3": SeqRecord(Seq(random_seq(10)), id="SEQ_3"),
    }

@pytest.fixture
def fasta_fn(tmpdir, sequences):
    fn = str(tmpdir / "sequences.fasta")
    SeqIO.write(sequences.values(), fn, "fasta")
    return fn

def write_metadata(tmpdir, metadata):
    fn = str(tmpdir / "metadata.tsv")
    with open(fn, "w") as fh:
        fh.write("\n".join(("\t".join(md) for md in metadata)))
    return fn

@pytest.fixture
def mock_priorities_file_valid(mocker):
    mocker.patch(
        "builtins.open", mocker.mock_open(read_data="strain1 5\nstrain2 6\nstrain3 8\n")
    )


@pytest.fixture
def mock_priorities_file_malformed(mocker):
    mocker.patch("builtins.open", mocker.mock_open(read_data="strain1 X\n"))


@pytest.fixture
def mock_run_shell_command(mocker):
    mocker.patch("augur.filter.run_shell_command")


@pytest.fixture
def mock_priorities_file_valid_with_spaces_and_tabs(mocker):
    mocker.patch(
        "builtins.open", mocker.mock_open(read_data="strain 1\t5\nstrain 2\t6\nstrain 3\t8\n")
    )

class TestFilter:
    def test_read_vcf_compressed(self):
        seq_keep, all_seq = augur.filter.read_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz"
        )

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_vcf_uncompressed(self):
        seq_keep, all_seq = augur.filter.read_vcf("tests/builds/tb/data/lee_2015.vcf")

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_priority_scores_valid(self, mock_priorities_file_valid):
        # builtins.open is stubbed, but we need a valid file to satisfy the existence check
        priorities = augur.filter.read_priority_scores(
            "tests/builds/tb/data/lee_2015.vcf"
        )

        assert priorities == {"strain1": 5, "strain2": 6, "strain3": 8}
        assert priorities["strain1"] == 5
        assert priorities["strain42"] == 0, "Default priority is 0 for unlisted sequences"

    def test_read_priority_scores_malformed(self, mock_priorities_file_malformed):
        with pytest.raises(ValueError):
            # builtins.open is stubbed, but we need a valid file to satisfy the existence check
            augur.filter.read_priority_scores("tests/builds/tb/data/lee_2015.vcf")

    def test_read_priority_scores_valid_with_spaces_and_tabs(self, mock_priorities_file_valid_with_spaces_and_tabs):
        # builtins.open is stubbed, but we need a valid file to satisfy the existence check
        priorities = augur.filter.read_priority_scores(
            "tests/builds/tb/data/lee_2015.vcf"
        )

        assert priorities == {"strain 1": 5, "strain 2": 6, "strain 3": 8}

    def test_read_priority_scores_does_not_exist(self):
        with pytest.raises(FileNotFoundError):
            augur.filter.read_priority_scores("/does/not/exist.txt")

    def test_write_vcf_compressed_input(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz", "output_file.vcf.gz", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --gzvcf tests/builds/tb/data/lee_2015.vcf.gz --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_input(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf.gz", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_compressed_output(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf.gz", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_output(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf", []
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )

    def test_write_vcf_dropped_samples(self, mock_run_shell_command):
        augur.filter.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf", ["x", "y", "z"]
        )

        augur.filter.run_shell_command.assert_called_once_with(
            "vcftools --remove-indv x --remove-indv y --remove-indv z --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )

    def test_filter_on_query_good(self, tmpdir, sequences):
        """Basic filter_on_query test"""
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        metadata, columns = read_metadata(meta_fn, as_data_frame=True)
        filtered = augur.filter.filter_by_query(set(sequences.keys()), metadata, 'quality=="good"')
        assert sorted(filtered) == ["SEQ_1", "SEQ_3"]

    def test_filter_on_query_subset(self, tmpdir):
        """Test filtering on query works when given fewer strains than metadata"""
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        metadata, columns = read_metadata(meta_fn, as_data_frame=True)
        filtered = augur.filter.filter_by_query({"SEQ_2"}, metadata, 'quality=="bad" & location=="colorado"')
        assert sorted(filtered) == ["SEQ_2"]

    def test_filter_run_with_query(self, tmpdir, fasta_fn, argparser):
        """Test that filter --query works as expected"""
        out_fn = str(tmpdir / "out.fasta")
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        args = argparser('-s %s --metadata %s -o %s --query "location==\'colorado\'"'
                         % (fasta_fn, meta_fn, out_fn))
        augur.filter.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_2"]

    def test_filter_run_with_query_and_include(self, tmpdir, fasta_fn, argparser):
        """Test that --include still works with filtering on query"""
        out_fn = str(tmpdir / "out.fasta")
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        include_fn = str(tmpdir / "include")
        open(include_fn, "w").write("SEQ_3")
        args = argparser('-s %s --metadata %s -o %s --query "quality==\'good\' & location==\'colorado\'" --include %s'
                         % (fasta_fn, meta_fn, out_fn, include_fn))
        augur.filter.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_3"]

    def test_filter_run_with_query_and_include_where(self, tmpdir, fasta_fn, argparser):
        """Test that --include_where still works with filtering on query"""
        out_fn = str(tmpdir / "out.fasta")
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        args = argparser('-s %s --metadata %s -o %s --query "quality==\'good\' & location==\'colorado\'" --include-where "location=nevada"'
                         % (fasta_fn, meta_fn, out_fn))
        augur.filter.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_3"]

    def test_filter_run_min_date(self, tmpdir, fasta_fn, argparser):
        """Test that filter --min-date is inclusive"""
        out_fn = str(tmpdir / "out.fasta")
        min_date = "2020-02-26"
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          ("SEQ_1","2020-02-XX"),
                                          ("SEQ_2","2020-02-26"),
                                          ("SEQ_3","2020-02-25")))
        args = argparser('-s %s --metadata %s -o %s --min-date %s'
                         % (fasta_fn, meta_fn, out_fn, min_date))
        augur.filter.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_2"]

    def test_filter_run_max_date(self, tmpdir, fasta_fn, argparser):
        """Test that filter --max-date is inclusive"""
        out_fn = str(tmpdir / "out.fasta")
        max_date = "2020-03-01"
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          ("SEQ_1","2020-03-XX"),
                                          ("SEQ_2","2020-03-01"),
                                          ("SEQ_3","2020-03-02")))
        args = argparser('-s %s --metadata %s -o %s --max-date %s'
                         % (fasta_fn, meta_fn, out_fn, max_date))
        augur.filter.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_2"]
