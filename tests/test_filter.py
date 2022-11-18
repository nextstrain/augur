import argparse
from textwrap import dedent
import numpy as np
import pandas as pd
import random
import shlex

import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from freezegun import freeze_time

import augur.filter
from augur.io.metadata import read_metadata
from augur.filter import FilterException

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
def mock_priorities_file_valid_with_spaces_and_tabs(mocker):
    mocker.patch(
        "builtins.open", mocker.mock_open(read_data="strain 1\t5\nstrain 2\t6\nstrain 3\t8\n")
    )

class TestFilter:
    def test_read_priority_scores_valid(self, mock_priorities_file_valid):
        # builtins.open is stubbed, but we need a valid file to satisfy the existence check
        priorities = augur.filter.read_priority_scores(
            "tests/builds/tb/data/lee_2015.vcf"
        )

        assert priorities == {"strain1": 5, "strain2": 6, "strain3": 8}
        assert priorities["strain1"] == 5
        assert priorities["strain42"] == -np.inf, "Default priority is negative infinity for unlisted sequences"

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

    def test_filter_on_query_good(self, tmpdir, sequences):
        """Basic filter_on_query test"""
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        metadata = read_metadata(meta_fn)
        filtered = augur.filter.filter_by_query(metadata, 'quality=="good"')
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
        augur.filter.run(args)
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
        augur.filter.run(args)
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
        augur.filter.run(args)
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
        assert e_info.value.__context__.message == dedent(f"""\
            Unable to determine date from '{argparse_value}'. Ensure it is in one of the supported formats:
            1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
            2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
            3. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
        """)


@pytest.fixture
def valid_metadata() -> pd.DataFrame:
    columns = ['strain', 'date', 'country']
    data = [
        ("SEQ_1","2020-01-XX","A"),
        ("SEQ_2","2020-02-01","A"),
        ("SEQ_3","2020-03-01","B"),
        ("SEQ_4","2020-04-01","B"),
        ("SEQ_5","2020-05-01","B")
    ]
    return pd.DataFrame.from_records(data, columns=columns).set_index('strain')

class TestFilterGroupBy:
    def test_filter_groupby_strain_subset(self, valid_metadata: pd.DataFrame):
        metadata = valid_metadata.copy()
        strains = ['SEQ_1', 'SEQ_3', 'SEQ_5']
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata)
        assert group_by_strain == {
            'SEQ_1': ('_dummy',),
            'SEQ_3': ('_dummy',),
            'SEQ_5': ('_dummy',)
        }
        assert skipped_strains == []

    def test_filter_groupby_dummy(self, valid_metadata: pd.DataFrame):
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata)
        assert group_by_strain == {
            'SEQ_1': ('_dummy',),
            'SEQ_2': ('_dummy',),
            'SEQ_3': ('_dummy',),
            'SEQ_4': ('_dummy',),
            'SEQ_5': ('_dummy',)
        }
        assert skipped_strains == []

    def test_filter_groupby_invalid_error(self, valid_metadata: pd.DataFrame):
        groups = ['invalid']
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['invalid']) were not found."

    def test_filter_groupby_invalid_warn(self, valid_metadata: pd.DataFrame, capsys):
        groups = ['country', 'year', 'month', 'invalid']
        metadata = valid_metadata.copy()
        strains = metadata.index.tolist()
        group_by_strain, _ = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1), 'unknown'),
            'SEQ_2': ('A', 2020, (2020, 2), 'unknown'),
            'SEQ_3': ('B', 2020, (2020, 3), 'unknown'),
            'SEQ_4': ('B', 2020, (2020, 4), 'unknown'),
            'SEQ_5': ('B', 2020, (2020, 5), 'unknown')
        }
        captured = capsys.readouterr()
        assert captured.err == "WARNING: Some of the specified group-by categories couldn't be found: invalid\nFiltering by group may behave differently than expected!\n"

    def test_filter_groupby_skip_ambiguous_year(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = "XXXX-02-01"
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_year', 'kwargs': ''}]

    def test_filter_groupby_skip_missing_date(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = None
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_year', 'kwargs': ''}]

    def test_filter_groupby_skip_ambiguous_month(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = "2020-XX-01"
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}]

    def test_filter_groupby_skip_missing_month(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata.at["SEQ_2", "date"] = "2020"
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 3)),
            'SEQ_4': ('B', 2020, (2020, 4)),
            'SEQ_5': ('B', 2020, (2020, 5))
        }
        assert skipped_strains == [{'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}]

    def test_filter_groupby_missing_year_error(self, valid_metadata: pd.DataFrame):
        groups = ['year']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['year']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_month_error(self, valid_metadata: pd.DataFrame):
        groups = ['month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['month']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_year_and_month_error(self, valid_metadata: pd.DataFrame):
        groups = ['year', 'month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        with pytest.raises(FilterException) as e_info:
            augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert str(e_info.value) == "The specified group-by categories (['year', 'month']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'."

    def test_filter_groupby_missing_date_warn(self, valid_metadata: pd.DataFrame, capsys):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata = metadata.drop('date', axis='columns')
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 'unknown', 'unknown'),
            'SEQ_2': ('A', 'unknown', 'unknown'),
            'SEQ_3': ('B', 'unknown', 'unknown'),
            'SEQ_4': ('B', 'unknown', 'unknown'),
            'SEQ_5': ('B', 'unknown', 'unknown')
        }
        captured = capsys.readouterr()
        assert captured.err == "WARNING: A 'date' column could not be found to group-by ['month', 'year'].\nFiltering by group may behave differently than expected!\n"
        assert skipped_strains == []

    def test_filter_groupby_no_strains(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        strains = []
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {}
        assert skipped_strains == []

    def test_filter_groupby_only_year_provided(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year']
        metadata = valid_metadata.copy()
        metadata['date'] = '2020'
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020),
            'SEQ_2': ('A', 2020),
            'SEQ_3': ('B', 2020),
            'SEQ_4': ('B', 2020),
            'SEQ_5': ('B', 2020)
        }
        assert skipped_strains == []

    def test_filter_groupby_month_with_only_year_provided(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata['date'] = '2020'
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {}
        assert skipped_strains == [
            {'strain': 'SEQ_1', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''},
            {'strain': 'SEQ_2', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''},
            {'strain': 'SEQ_3', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''},
            {'strain': 'SEQ_4', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''},
            {'strain': 'SEQ_5', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}
        ]

    def test_filter_groupby_only_year_month_provided(self, valid_metadata: pd.DataFrame):
        groups = ['country', 'year', 'month']
        metadata = valid_metadata.copy()
        metadata['date'] = '2020-01'
        strains = metadata.index.tolist()
        group_by_strain, skipped_strains = augur.filter.get_groups_for_subsampling(strains, metadata, group_by=groups)
        assert group_by_strain == {
            'SEQ_1': ('A', 2020, (2020, 1)),
            'SEQ_2': ('A', 2020, (2020, 1)),
            'SEQ_3': ('B', 2020, (2020, 1)),
            'SEQ_4': ('B', 2020, (2020, 1)),
            'SEQ_5': ('B', 2020, (2020, 1))
        }
        assert skipped_strains == []
