import argparse
import pytest
import shlex
import sqlite3
from freezegun import freeze_time
from textwrap import dedent

import augur.filter
from augur.dates import (
    any_to_numeric_type_min,
    any_to_numeric_type_max,
)
from augur.filter_support.exceptions import FilterException
from augur.filter_support.subsample import get_valid_group_by_cols
from augur.filter_support.db.sqlite import (
    FilterSQLite,
    EXCLUDE_COL,
    INCLUDE_COL,
    FILTER_REASON_COL,
    FILTER_REASON_KWARGS_COL,
    GROUP_SIZES_TABLE_NAME,
    METADATA_FILTER_REASON_TABLE_NAME,
    METADATA_TABLE_NAME,
    OUTPUT_METADATA_TABLE_NAME,
    PRIORITIES_TABLE_NAME,
    SEQUENCE_INDEX_TABLE_NAME,
)


def parse_args(args:str):
    parser = argparse.ArgumentParser()
    augur.filter.register_arguments(parser)
    return parser.parse_args(shlex.split(args))


def write_file(tmpdir, filename:str, content:str):
    filepath = str(tmpdir / filename)
    with open(filepath, "w") as handle:
        handle.write(content)
    return filepath


def write_metadata(tmpdir, metadata):
    content = "\n".join(("\t".join(md) for md in metadata))
    return write_file(tmpdir, "metadata.tsv", content)


def get_filter_obj_run(args:argparse.Namespace):
    """Returns a filter object connected to an in-memory database with run() invoked."""
    # use an in-memory database for tests since:
    # 1. test data is not large
    # 2. in-memory I/O is generally faster
    obj = FilterSQLite(args, in_memory_db=True)
    obj.run()
    return obj


def get_valid_args(data, tmpdir):
    """Returns an argparse.Namespace with metadata and output_strains"""
    meta_fn = write_metadata(tmpdir, data)
    return parse_args(f'--metadata {meta_fn} --output-strains {tmpdir / "strains.txt"}')


def query_fetchall(filter_obj:FilterSQLite, query:str):
    with filter_obj.get_db_context() as con:
        return con.execute(query).fetchall()


def query_fetchall_dict(filter_obj:FilterSQLite, query:str):
    with filter_obj.get_db_context() as con:
        con.row_factory = sqlite3.Row
        return [dict(row) for row in con.execute(query)]


class TestFiltering:
    def test_filter_by_query(self, tmpdir):
        """Filter by a query expresssion."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.query = 'quality=="good"'
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_query'
        """)
        assert results == [("SEQ_2",)]

    def test_filter_by_query_two_conditions(self, tmpdir):
        """Filter by a query expresssion with two conditions."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.query = 'quality=="good" AND location=="colorado"'
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_query'
        """)
        assert results == [("SEQ_2",), ("SEQ_3",)]

    def test_filter_by_query_and_include_strains(self, tmpdir):
        """Filter by a query expresssion and force-include a strain."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        include_fn = str(tmpdir / "include.txt")
        open(include_fn, "w").write("SEQ_3")
        args = get_valid_args(data, tmpdir)
        args.query = 'quality=="good" AND location=="colorado"'
        args.include = [include_fn]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE NOT {EXCLUDE_COL} OR {INCLUDE_COL}
        """)
        assert results == [("SEQ_1",), ("SEQ_3",)]

    def test_filter_by_query_and_include_where(self, tmpdir):
        """Filter by a query expresssion and force-include a strain."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.query = 'quality=="good" AND location=="colorado"'
        args.include_where = ['location=nevada']
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE NOT {EXCLUDE_COL} OR {INCLUDE_COL}
        """)
        assert results == [("SEQ_1",), ("SEQ_3",)]

    def test_filter_by_min_date(self, tmpdir):
        """Filter by min date, inclusive."""
        data = [("strain","date"),
                ("SEQ_1","2020-02-XX"),
                ("SEQ_2","2020-02-26"),
                ("SEQ_3","2020-02-25"),
                ("SEQ_4","?")]
        args = get_valid_args(data, tmpdir)
        args.min_date = any_to_numeric_type_min('2020-02-26')
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_min_date'
        """)
        assert results == [("SEQ_3",), ("SEQ_4",)]

    def test_filter_by_max_date(self, tmpdir):
        """Filter by max date, inclusive."""
        data = [("strain","date"),
                ("SEQ_1","2020-03-XX"),
                ("SEQ_2","2020-03-01"),
                ("SEQ_3","2020-03-02"),
                ("SEQ_4","?")]
        args = get_valid_args(data, tmpdir)
        args.max_date = any_to_numeric_type_max('2020-03-01')
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_max_date'
        """)
        assert results == [("SEQ_3",), ("SEQ_4",)]

    def test_filter_incomplete_year(self, tmpdir):
        """Test that 2020 is evaluated as 2020-XX-XX"""
        data = [("strain","date"),
                ("SEQ_1","2020.0"),
                ("SEQ_2","2020"),
                ("SEQ_3","2020-XX-XX")]
        args = get_valid_args(data, tmpdir)
        args.min_date = any_to_numeric_type_min('2020-02-01')
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_min_date'
        """)
        assert results == [("SEQ_1",)]

    def test_filter_date_formats(self, tmpdir):
        """Test that 2020.0, 2020, and 2020-XX-XX all pass --min-date 2019"""
        data = [("strain","date"),
                ("SEQ_1","2020.0"),
                ("SEQ_2","2020"),
                ("SEQ_3","2020-XX-XX")]
        args = get_valid_args(data, tmpdir)
        args.min_date = any_to_numeric_type_min('2019')
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_min_date'
        """)
        assert results == []

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
    def test_filter_relative_dates(self, tmpdir, argparse_params, metadata_rows, output_sorted_expected):
        """Test that various relative dates work"""
        out_fn = str(tmpdir / "filtered.txt")
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          *metadata_rows))
        args = parse_args(f'--metadata {meta_fn} --output-strains {out_fn} {argparse_params}')
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
    def test_filter_relative_dates_error(self, tmpdir, argparse_flag, argparse_value):
        """Test that invalid dates fail"""
        out_fn = str(tmpdir / "filtered.txt")
        meta_fn = write_metadata(tmpdir, (("strain","date"),
                                          ("SEQ_1","2020-03-23")))
        with pytest.raises(SystemExit) as e_info:
            parse_args(f'--metadata {meta_fn} --output-strains {out_fn} {argparse_flag} {argparse_value}')
        assert e_info.value.__context__.message == dedent(f"""\
            Unable to determine date from '{argparse_value}'. Ensure it is in one of the supported formats:
            1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
            2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
            3. an ambiguous date in ISO 8601-like format (e.g. '2020-06-XX', '2020-XX-XX') or
            4. an incomplete date in ISO 8601-like format (e.g. '2020-06', '2020') or
            5. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
        """)

    def test_filter_by_ambiguous_date_year(self, tmpdir):
        """Filter out dates with ambiguous year."""
        data = [("strain","date"),
                ("SEQ_1","XXXX"),
                ("SEQ_2","2020-XX"),
                ("SEQ_3","2020-03-XX"),
                ("SEQ_4","?")]
        args = get_valid_args(data, tmpdir)
        args.exclude_ambiguous_dates_by = 'year'
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain, {FILTER_REASON_KWARGS_COL}
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_ambiguous_date'
        """)
        assert results == [
            ('SEQ_1', '[["ambiguity", "year"]]'),
            ('SEQ_4', '[["ambiguity", "year"]]')
        ]

    def test_filter_by_ambiguous_date_month(self, tmpdir):
        """Filter out dates with ambiguous month."""
        data = [("strain","date"),
                ("SEQ_1","XXXX"),
                ("SEQ_2","2020-XX"),
                ("SEQ_3","2020-03-XX"),
                ("SEQ_4","?")]
        args = get_valid_args(data, tmpdir)
        args.exclude_ambiguous_dates_by = 'month'
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain, {FILTER_REASON_KWARGS_COL}
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_ambiguous_date'
        """)
        assert results == [
            ('SEQ_1', '[["ambiguity", "month"]]'),
            ('SEQ_2', '[["ambiguity", "month"]]'),
            ('SEQ_4', '[["ambiguity", "month"]]')
        ]

    def test_filter_by_ambiguous_date_day(self, tmpdir):
        """Filter out dates with ambiguous day."""
        data = [("strain","date"),
                ("SEQ_1","XXXX"),
                ("SEQ_2","2020-XX"),
                ("SEQ_3","2020-03-XX"),
                ("SEQ_4","2020-03-02"),
                ("SEQ_5","?")]
        args = get_valid_args(data, tmpdir)
        args.exclude_ambiguous_dates_by = 'day'
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain, {FILTER_REASON_KWARGS_COL}
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_ambiguous_date'
        """)
        assert results == [
            ('SEQ_1', '[["ambiguity", "day"]]'),
            ('SEQ_2', '[["ambiguity", "day"]]'),
            ('SEQ_3', '[["ambiguity", "day"]]'),
            ('SEQ_5', '[["ambiguity", "day"]]')
        ]

    def test_filter_by_ambiguous_date_any(self, tmpdir):
        """Filter out dates with any ambiguity."""
        data = [("strain","date"),
                ("SEQ_1","XXXX"),
                ("SEQ_2","2020-XX"),
                ("SEQ_3","2020-03-XX"),
                ("SEQ_4","2020-03-02"),
                ("SEQ_5","?")]
        args = get_valid_args(data, tmpdir)
        args.exclude_ambiguous_dates_by = 'any'
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain, {FILTER_REASON_KWARGS_COL}
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_ambiguous_date'
        """)
        assert results == [
            ('SEQ_1', '[["ambiguity", "any"]]'),
            ('SEQ_2', '[["ambiguity", "any"]]'),
            ('SEQ_3', '[["ambiguity", "any"]]'),
            ('SEQ_5', '[["ambiguity", "any"]]')
        ]

    def test_filter_by_exclude_where(self, tmpdir):
        """Filter by an expression that matches location equal to colorado."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.exclude_where = ["location=colorado"]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_exclude_where'
        """)
        assert results == [("SEQ_1",), ("SEQ_2",)]

    def test_filter_by_exclude_where_missing_column_error(self, tmpdir):
        """Try filtering by an expression matching on an invalid column."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.exclude_where = ["invalid=colorado"]
        with pytest.raises(FilterException) as e_info:
            get_filter_obj_run(args)
        assert str(e_info.value) == 'no such column: metadata.invalid'

    def test_force_include_where_missing_column_error(self, tmpdir):
        """Try filtering by an expression matching on an invalid column."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.include_where = ["invalid=colorado"]
        with pytest.raises(FilterException) as e_info:
            get_filter_obj_run(args)
        assert str(e_info.value) == 'no such column: metadata.invalid'

    def test_filter_by_min_length(self, tmpdir):
        """Filter by minimum sequence length of 3."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",)]
        args = get_valid_args(data, tmpdir)
        fasta_lines = [
            ">SEQ_1", "aa",
            ">SEQ_2", "aaa",
            ">SEQ_3", "nnnn",
        ]
        args.sequences = write_file(tmpdir, "sequences.fasta", "\n".join(fasta_lines))
        args.min_length = 3
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_sequence_length'
        """)
        assert results == [("SEQ_1",), ("SEQ_3",)]

    def test_filter_by_non_nucleotide(self, tmpdir):
        """Filter out sequences with at least 1 invalid nucleotide character."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",),
                ("SEQ_4",)]
        args = get_valid_args(data, tmpdir)
        fasta_lines = [
            ">SEQ_1", "aaaa",
            ">SEQ_2", "nnnn",
            ">SEQ_3", "xxxx",
            ">SEQ_4", "aaax",
        ]
        args.sequences = write_file(tmpdir, "sequences.fasta", "\n".join(fasta_lines))
        args.non_nucleotide = True
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_non_nucleotide'
        """)
        assert results == [("SEQ_3",), ("SEQ_4",)]

    def test_filter_by_exclude_all(self, tmpdir):
        """Filter out all sequences."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",),
                ("SEQ_4",)]
        args = get_valid_args(data, tmpdir)
        args.exclude_all = True
        with pytest.raises(FilterException) as e_info:
            get_filter_obj_run(args)
        assert str(e_info.value) == "All samples have been dropped! Check filter rules and metadata file format."

    def test_filter_by_exclude_strains(self, tmpdir):
        """Exclude strains from a file."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",),
                ("SEQ_4",)]
        args = get_valid_args(data, tmpdir)
        exclude_seqs = [
            "SEQ_1",
            "SEQ_3",
        ]
        args.exclude = [write_file(tmpdir, "exclude.txt", "\n".join(exclude_seqs))]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_exclude_strains'
        """)
        assert results == [("SEQ_1",), ("SEQ_3",)]

    def test_filter_by_force_include_strains(self, tmpdir):
        """Force-include strains from a file."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",),
                ("SEQ_4",)]
        args = get_valid_args(data, tmpdir)
        include_seqs = [
            "SEQ_1",
            "SEQ_3",
        ]
        args.include = [write_file(tmpdir, "include.txt", "\n".join(include_seqs))]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'force_include_strains'
        """)
        assert results == [("SEQ_1",), ("SEQ_3",)]

    def test_filter_by_exclude_where_negative(self, tmpdir):
        """Filter by an expression that matches location not equal to colorado."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        args.exclude_where = ["location!=colorado"]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain, {FILTER_REASON_KWARGS_COL}
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'filter_by_exclude_where'
        """)
        assert results == [('SEQ_3', '[["exclude_where", "location!=colorado"]]')]


class TestFilterGroupBy:
    def test_filter_groupby_invalid_error(self):
        metadata_cols = {"strain", "date", "country"}
        groups = ["invalid"]
        with pytest.raises(FilterException) as e_info:
            get_valid_group_by_cols(groups, metadata_cols)
        assert str(e_info.value) == "The specified group-by categories (['invalid']) were not found."

    def test_filter_groupby_invalid_warn(self, capsys):
        metadata_cols = {"strain", "date", "country"}
        groups = ["country", "year", "month", "invalid"]
        valid_group_by_cols = get_valid_group_by_cols(groups, metadata_cols)
        assert valid_group_by_cols == ["country", "year", "month"]
        captured = capsys.readouterr()
        assert captured.err == dedent("""\
            WARNING: Some of the specified group-by categories couldn't be found: invalid
            Filtering by group may behave differently than expected!
        """)

    def test_filter_groupby_skip_ambiguous_year(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01-XX","A"),
            ("SEQ_2","XXXX","A"),
            ("SEQ_3","2020-03-01","B"),
            ("SEQ_4","2020-04-01","B"),
            ("SEQ_5","2020-05-01","B")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'skip_group_by_with_ambiguous_year'
        """)
        assert results == [("SEQ_2",)]

    def test_filter_groupby_skip_missing_date(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01-XX","A"),
            ("SEQ_2","","A"),
            ("SEQ_3","2020-03-01","B"),
            ("SEQ_4","2020-04-01","B"),
            ("SEQ_5","2020-05-01","B")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year", "month"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'skip_group_by_with_ambiguous_month'
        """)
        assert results == [("SEQ_2",)]

    def test_filter_groupby_skip_ambiguous_month(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01-XX","A"),
            ("SEQ_2","2020-XX-XX","A"),
            ("SEQ_3","2020-03-01","B"),
            ("SEQ_4","2020-04-01","B"),
            ("SEQ_5","2020-05-01","B")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year", "month"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'skip_group_by_with_ambiguous_month'
        """)
        assert results == [("SEQ_2",)]

    def test_filter_groupby_skip_missing_month(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01-XX","A"),
            ("SEQ_2","2020","A"),
            ("SEQ_3","2020-03-01","B"),
            ("SEQ_4","2020-04-01","B"),
            ("SEQ_5","2020-05-01","B")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year", "month"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'skip_group_by_with_ambiguous_month'
        """)
        assert results == [("SEQ_2",)]

    def test_filter_groupby_missing_year_error(self):
        metadata_cols = {"strain", "country"}
        groups = ["year"]
        with pytest.raises(FilterException) as e_info:
            get_valid_group_by_cols(groups, metadata_cols)
        assert str(e_info.value) == "The specified group-by categories (['year']) were not found. Note that using 'year' or 'year month' requires a column called 'date'."

    def test_filter_groupby_missing_month_error(self):
        metadata_cols = {"strain", "country"}
        groups = ["month"]
        with pytest.raises(FilterException) as e_info:
            get_valid_group_by_cols(groups, metadata_cols)
        assert str(e_info.value) == "The specified group-by categories (['month']) were not found. Note that using 'year' or 'year month' requires a column called 'date'."

    def test_filter_groupby_missing_year_and_month_error(self):
        metadata_cols = {"strain", "country"}
        groups = ["year", "month"]
        with pytest.raises(FilterException) as e_info:
            get_valid_group_by_cols(groups, metadata_cols)
        assert str(e_info.value) == "The specified group-by categories (['year', 'month']) were not found. Note that using 'year' or 'year month' requires a column called 'date'."

    def test_filter_groupby_missing_date_warn(self, capsys):
        metadata_cols = {"strain", "country"}
        groups = ["country", "year", "month"]
        valid_group_by_cols = get_valid_group_by_cols(groups, metadata_cols)
        assert valid_group_by_cols == ["country"]
        captured = capsys.readouterr()
        assert captured.err == dedent("""\
            WARNING: A 'date' column could not be found to group-by year.
            WARNING: A 'date' column could not be found to group-by month.
            WARNING: Some of the specified group-by categories couldn't be found: year, month
            Filtering by group may behave differently than expected!
        """)

    def test_filter_groupby_only_year_provided(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01-XX","A"),
            ("SEQ_2","2020","B"),
            ("SEQ_3","2020-03-01","C"),
            ("SEQ_4","2020-04-01","C"),
            ("SEQ_5","2020-05-01","C")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT country, year FROM {GROUP_SIZES_TABLE_NAME}
        """)
        assert results == [
            ("A", 2020),
            ("B", 2020),
            ("C", 2020)
        ]

    def test_filter_groupby_month_with_only_year_provided(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01-XX","A"),
            ("SEQ_2","2020","B"),
            ("SEQ_3","2020-03-01","C"),
            ("SEQ_4","2020-04-01","C"),
            ("SEQ_5","2020-05-01","C")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year", "month"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT country, year, month FROM {GROUP_SIZES_TABLE_NAME}
        """)
        assert results == [
            ("A", 2020, 1),
            ("C", 2020, 3),
            ("C", 2020, 4),
            ("C", 2020, 5)
        ]
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} = 'skip_group_by_with_ambiguous_month'
        """)
        assert results == [("SEQ_2",)]

    def test_filter_groupby_only_year_month_provided(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020-01","A"),
            ("SEQ_2","2020-02","B"),
            ("SEQ_3","2020-03","C"),
            ("SEQ_4","2020-04","D"),
            ("SEQ_5","2020-05","E")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year", "month"]
        args.sequences_per_group = 1
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT country, year, month FROM {GROUP_SIZES_TABLE_NAME}
        """)
        assert results == [
            ("A", 2020, 1),
            ("B", 2020, 2),
            ("C", 2020, 3),
            ("D", 2020, 4),
            ("E", 2020, 5)
        ]
        results = query_fetchall(filter_obj, f"""
            SELECT strain FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} IS NULL
        """)
        assert results == [("SEQ_1",), ("SEQ_2",), ("SEQ_3",), ("SEQ_4",), ("SEQ_5",)]

    def test_all_samples_dropped(self, tmpdir):
        data = [
            ("strain","date","country"),
            ("SEQ_1","2020","A"),
            ("SEQ_2","2020","B"),
            ("SEQ_3","2020","C"),
            ("SEQ_4","2020","D"),
            ("SEQ_5","2020","E")
        ]
        args = get_valid_args(data, tmpdir)
        args.group_by = ["country", "year", "month"]
        args.sequences_per_group = 1
        with pytest.raises(FilterException) as e_info:
            get_filter_obj_run(args)
        assert str(e_info.value) == "All samples have been dropped! Check filter rules and metadata file format."


class TestDataLoading:
    def test_load_metadata(self, tmpdir):
        """Load a metadata file."""
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good"),
                ("SEQ_2","colorado","bad"),
                ("SEQ_3","nevada","good")]
        args = get_valid_args(data, tmpdir)
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"SELECT * FROM {METADATA_TABLE_NAME}")
        assert [row[1:] for row in results] == data[1:]

    def test_load_priority_scores_valid(self, tmpdir):
        """Load a priority score file."""
        content = "strain1\t5\nstrain2\t6\nstrain3\t8\n"
        filter_obj = get_filter_obj_with_priority_loaded(tmpdir, content)
        filter_obj.db_load_priorities_table()
        results = query_fetchall(filter_obj, f"SELECT * FROM {PRIORITIES_TABLE_NAME}")
        assert results == [(0, "strain1", 5.0), (1, "strain2", 6.0), (2, "strain3", 8.0)]

    @pytest.mark.skip(reason="this isn't trivial with SQLite's flexible typing rules")
    def test_load_priority_scores_malformed(self, tmpdir):
        """Attempt to load a priority score file with non-float in priority column raises a ValueError."""
        content = "strain1 X\n"
        filter_obj = get_filter_obj_with_priority_loaded(tmpdir, content)
        with pytest.raises(ValueError) as e_info:
            filter_obj.db_load_priorities_table()
        assert str(e_info.value) == f"Failed to parse priority file {filter_obj.args.priority}."

    def test_load_priority_scores_valid_with_spaces_and_tabs(self, tmpdir):
        """Load a priority score file with spaces in strain names."""
        content = "strain 1\t5\nstrain 2\t6\nstrain 3\t8\n"
        filter_obj = get_filter_obj_with_priority_loaded(tmpdir, content)
        filter_obj.db_load_priorities_table()
        results = query_fetchall(filter_obj, f"SELECT * FROM {PRIORITIES_TABLE_NAME}")
        assert results == [(0, "strain 1", 5.0), (1, "strain 2", 6.0), (2, "strain 3", 8.0)]

    def test_load_priority_scores_does_not_exist(self, tmpdir):
        """Attempt to load a non-existant priority score file raises a FileNotFoundError."""
        invalid_priorities_fn = str(tmpdir / "does/not/exist.txt")
        # metadata is a required arg but we don't need it
        data = [("strain","location","quality"),
                ("SEQ_1","colorado","good")]
        args = get_valid_args(data, tmpdir)
        args.priority = invalid_priorities_fn
        filter_obj = get_filter_obj_run(args)
        with pytest.raises(FileNotFoundError):
            filter_obj.db_load_priorities_table()

    def test_load_invalid_id_column(self, tmpdir):
        data = [
            ("invalid_name","date","country"),
            ("SEQ_1","2020-01-XX","A"),
        ]
        args = get_valid_args(data, tmpdir)
        with pytest.raises(ValueError) as e_info:
            get_filter_obj_run(args)
        assert str(e_info.value) == "None of the possible id columns (['strain', 'name']) were found in the metadata's columns ('invalid_name', 'date', 'country')"

    def test_load_custom_id_column(self, tmpdir):
        data = [
            ("custom_id_col","date","country"),
            ("SEQ_1","2020-01-XX","A"),
        ]
        args = get_valid_args(data, tmpdir)
        args.metadata_id_columns = ["custom_id_col"]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT custom_id_col FROM {METADATA_TABLE_NAME}
        """)
        assert results == [("SEQ_1",)]

    def test_load_custom_id_column_with_spaces(self, tmpdir):
        data = [
            ("strain name with spaces","date","country"),
            ("SEQ_1","2020-01-XX","A"),
        ]
        args = get_valid_args(data, tmpdir)
        args.metadata_id_columns = ["strain name with spaces"]
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"""
            SELECT "strain name with spaces" FROM {METADATA_TABLE_NAME}
        """)
        assert results == [("SEQ_1",)]

    def test_load_priority_scores_extra_column(self, tmpdir):
        """Attempt to load a priority score file with an extra column raises a ValueError."""
        content = "strain1\t5\tbad_col\n"
        filter_obj = get_filter_obj_with_priority_loaded(tmpdir, content)
        with pytest.raises(ValueError) as e_info:
            filter_obj.db_load_priorities_table()
        assert str(e_info.value) == f"Failed to parse priority file {filter_obj.args.priority}."

    def test_load_priority_scores_missing_column(self, tmpdir):
        """Attempt to load a priority score file with a missing column raises a ValueError."""
        content = "strain1\n"
        filter_obj = get_filter_obj_with_priority_loaded(tmpdir, content)
        with pytest.raises(ValueError) as e_info:
            filter_obj.db_load_priorities_table()
        assert str(e_info.value) == f"Failed to parse priority file {filter_obj.args.priority}."

    def test_load_sequences_subset_strains(self, tmpdir):
        """Loading sequences filters output to the intersection of strains from metadata and sequences."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",)]
        args = get_valid_args(data, tmpdir)
        fasta_lines = [
            ">SEQ_1", "aaaa",
            ">SEQ_3", "aaaa",
            ">SEQ_4", "nnnn",
        ]
        args.sequences = write_file(tmpdir, "sequences.fasta", "\n".join(fasta_lines))
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"SELECT strain FROM {OUTPUT_METADATA_TABLE_NAME}")
        assert results == [("SEQ_1",), ("SEQ_3",)]

    def test_generate_sequence_index(self, tmpdir):
        """Loading sequences filters output to the intersection of strains from metadata and sequences."""
        data = [("strain",),
                ("SEQ_1",),
                ("SEQ_2",),
                ("SEQ_3",)]
        args = get_valid_args(data, tmpdir)
        fasta_lines = [
            ">SEQ_1", "aaaa",
            ">SEQ_3", "aaaa",
            ">SEQ_4", "nnnn",
        ]
        args.sequences = write_file(tmpdir, "sequences.fasta", "\n".join(fasta_lines))
        filter_obj = get_filter_obj_run(args)
        results = query_fetchall(filter_obj, f"SELECT strain, A, N FROM {SEQUENCE_INDEX_TABLE_NAME}")
        print(results)
        assert results == [("SEQ_1", 4, 0), ("SEQ_3", 4, 0), ("SEQ_4", 0, 4)]


def get_filter_obj_with_priority_loaded(tmpdir, content:str):
    priorities_fn = write_file(tmpdir, "priorities.txt", content)
    # metadata is a required arg but we don't need it
    data = [("strain","location","quality"),
            ("SEQ_1","colorado","good")]
    args = get_valid_args(data, tmpdir)
    args.priority = priorities_fn
    return get_filter_obj_run(args)
