import argparse
import random
import shlex

import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from freezegun import freeze_time

import augur.filter
import augur.filter._run
import augur.filter.include_exclude_rules
from augur.io.metadata import read_metadata

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


class TestFilter:
    def test_filter_on_query_good(self, tmpdir, sequences):
        """Basic filter_on_query test"""
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        metadata = read_metadata(meta_fn)
        filtered = augur.filter.include_exclude_rules.filter_by_query(metadata, 'quality=="good"')
        assert sorted(filtered) == ["SEQ_1", "SEQ_3"]

    def test_filter_run_with_query(self, tmpdir, fasta_fn, argparser):
        """Test that filter --query works as expected"""
        out_fn = str(tmpdir / "out.fasta")
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        args = argparser('-s %s --metadata %s -o %s --query "location==\'colorado\'"'
                         % (fasta_fn, meta_fn, out_fn))
        augur.filter._run.run(args)
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
        augur.filter._run.run(args)
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
        augur.filter._run.run(args)
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
        augur.filter._run.run(args)
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
        augur.filter._run.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_2"]

    def test_filter_incomplete_year(self, tmpdir, fasta_fn, argparser):
        """Test that 2020 is evaluated as 2020-XX-XX"""
        out_fn = str(tmpdir / "out.fasta")
        min_date = "2020-02-01"
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          ("SEQ_1","2020.0"),
                                          ("SEQ_2","2020"),
                                          ("SEQ_3","2020-XX-XX")))
        args = argparser('-s %s --metadata %s -o %s --min-date %s'
                         % (fasta_fn, meta_fn, out_fn, min_date))
        augur.filter._run.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_2", "SEQ_3"]

    def test_filter_date_formats(self, tmpdir, fasta_fn, argparser):
        """Test that 2020.0, 2020, and 2020-XX-XX all pass --min-date 2019"""
        out_fn = str(tmpdir / "out.fasta")
        min_date = "2019"
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          ("SEQ_1","2020.0"),
                                          ("SEQ_2","2020"),
                                          ("SEQ_3","2020-XX-XX")))
        args = argparser('-s %s --metadata %s -o %s --min-date %s'
                         % (fasta_fn, meta_fn, out_fn, min_date))
        augur.filter._run.run(args)
        output = SeqIO.to_dict(SeqIO.parse(out_fn, "fasta"))
        assert list(output.keys()) == ["SEQ_1", "SEQ_2", "SEQ_3"]

    @freeze_time("2020-03-25")
    @pytest.mark.parametrize(
        "argparse_params, metadata_rows, output_sorted_expected",
        [
            (
                "--min-date 1D",
                (
                    ("SEQ_1","2020-03-23"),
                    ("SEQ_2","2020-03-24"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_2", "SEQ_3"],
            ),
            (
                "--max-date 1D",
                (
                    ("SEQ_1","2020-03-23"),
                    ("SEQ_2","2020-03-24"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_1", "SEQ_2"],
            ),
            (
                "--min-date 4W",
                (
                    ("SEQ_1","2020-02-25"),
                    ("SEQ_2","2020-02-26"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_2", "SEQ_3"],
            ),
            (
                "--max-date 4W",
                (
                    ("SEQ_1","2020-02-25"),
                    ("SEQ_2","2020-02-26"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_1", "SEQ_2"],
            ),
            (
                "--min-date 1M",
                (
                    ("SEQ_1","2020-01-25"),
                    ("SEQ_2","2020-02-25"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_2", "SEQ_3"],
            ),
            (
                "--max-date 1M",
                (
                    ("SEQ_1","2020-01-25"),
                    ("SEQ_2","2020-02-25"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_1", "SEQ_2"],
            ),
            (
                "--min-date P1M",
                (
                    ("SEQ_1","2020-01-25"),
                    ("SEQ_2","2020-02-25"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_2", "SEQ_3"],
            ),
            (
                "--max-date P1M",
                (
                    ("SEQ_1","2020-01-25"),
                    ("SEQ_2","2020-02-25"),
                    ("SEQ_3","2020-03-25"),
                ),
                ["SEQ_1", "SEQ_2"],
            ),
            (
                "--min-date 2Y",
                (
                    ("SEQ_1","2017-03-25"),
                    ("SEQ_2","2018-03-25"),
                    ("SEQ_3","2019-03-25"),
                ),
                ["SEQ_2", "SEQ_3"],
            ),
            (
                "--max-date 2Y",
                (
                    ("SEQ_1","2017-03-25"),
                    ("SEQ_2","2018-03-25"),
                    ("SEQ_3","2019-03-25"),
                ),
                ["SEQ_1", "SEQ_2"],
            ),
            (
                "--min-date 1Y2W5D",
                (
                    ("SEQ_1","2019-03-05"),
                    ("SEQ_2","2019-03-06"),
                    ("SEQ_3","2019-03-07"),
                ),
                ["SEQ_2", "SEQ_3"],
            ),
            (
                "--max-date 1Y2W5D",
                (
                    ("SEQ_1","2019-03-05"),
                    ("SEQ_2","2019-03-06"),
                    ("SEQ_3","2019-03-07"),
                ),
                ["SEQ_1", "SEQ_2"],
            ),
        ],
    )
    def test_filter_relative_dates(self, tmpdir, argparser, argparse_params, metadata_rows, output_sorted_expected):
        """Test that various relative dates work"""
        out_fn = str(tmpdir / "filtered.txt")
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          *metadata_rows))
        args = argparser(f'--metadata {meta_fn} --output-strains {out_fn} {argparse_params}')
        augur.filter._run.run(args)
        with open(out_fn) as f:
            output_sorted = sorted(line.rstrip() for line in f)
        assert output_sorted == output_sorted_expected

    @freeze_time("2020-03-25")
    @pytest.mark.parametrize(
        "argparse_flag, argparse_value",
        [
            ("--min-date", "3000Y"),
            ("--max-date", "3000Y"),
            ("--min-date", "invalid"),
            ("--max-date", "invalid"),
        ],
    )
    def test_filter_relative_dates_error(self, tmpdir, argparser, argparse_flag, argparse_value):
        """Test that invalid dates fail"""
        out_fn = str(tmpdir / "filtered.txt")
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          ("SEQ_1","2020-03-23")))
        with pytest.raises(SystemExit) as e_info:
            argparser(f'--metadata {meta_fn} --output-strains {out_fn} {argparse_flag} {argparse_value}')
        assert f"Invalid date '{argparse_value}'" in e_info.value.__context__.message
