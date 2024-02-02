import json
from augur.errors import AugurError
from augur.io.print import print_err
from augur.io.sqlite3 import Sqlite3Database
from augur.types import EmptyOutputReportingMethod
from . import constants, include_exclude_rules


def print_report(args):
    """Print a report of how many strains were dropped and reasoning."""
    # Calculate the number of strains that don't exist in either metadata or
    # sequences.
    num_excluded_by_lack_of_metadata = _get_num_excluded_by_lack_of_metadata()

    # Calculate the number of strains passed and filtered.
    total_strains_passed = _get_total_strains_passed()
    total_strains_filtered = _get_num_metadata_strains() + num_excluded_by_lack_of_metadata - total_strains_passed

    print_err(f"{total_strains_filtered} {'strain was' if total_strains_filtered == 1 else 'strains were'} dropped during filtering")

    if num_excluded_by_lack_of_metadata:
        print_err(f"\t{num_excluded_by_lack_of_metadata} had no metadata")

    report_template_by_filter_name = {
        include_exclude_rules.filter_by_sequence_index.__name__: "{count} had no sequence data",
        include_exclude_rules.filter_by_exclude_all.__name__: "{count} {were} dropped by `--exclude-all`",
        include_exclude_rules.filter_by_exclude.__name__: "{count} {were} dropped because {they} {were} in {exclude_file}",
        include_exclude_rules.filter_by_exclude_where.__name__: "{count} {were} dropped because of '{exclude_where}'",
        include_exclude_rules.filter_by_query.__name__: "{count} {were} filtered out by the query: \"{query}\"",
        include_exclude_rules.filter_by_ambiguous_date.__name__: "{count} {were} dropped because of their ambiguous date in {ambiguity}",
        include_exclude_rules.filter_by_min_date.__name__: "{count} {were} dropped because {they} {were} earlier than {min_date} or missing a date",
        include_exclude_rules.filter_by_max_date.__name__: "{count} {were} dropped because {they} {were} later than {max_date} or missing a date",
        include_exclude_rules.filter_by_min_length.__name__: "{count} {were} dropped because {they} {were} shorter than the minimum length of {min_length}bp when only counting standard nucleotide characters A, C, G, or T (case-insensitive)",
        include_exclude_rules.filter_by_max_length.__name__: "{count} {were} dropped because {they} {were} longer than the maximum length of {max_length}bp when only counting standard nucleotide characters A, C, G, or T (case-insensitive)",
        include_exclude_rules.filter_by_non_nucleotide.__name__: "{count} {were} dropped because {they} had non-nucleotide characters",
        include_exclude_rules.skip_group_by_with_ambiguous_year.__name__: "{count} {were} dropped during grouping due to ambiguous year information",
        include_exclude_rules.skip_group_by_with_ambiguous_month.__name__: "{count} {were} dropped during grouping due to ambiguous month information",
        include_exclude_rules.skip_group_by_with_ambiguous_day.__name__: "{count} {were} dropped during grouping due to ambiguous day information",
        include_exclude_rules.force_include_strains.__name__: "{count} {were} force-included because {they} {were} in {include_file}",
        include_exclude_rules.force_include_where.__name__: "{count} {were} force-included because of '{include_where}'",
    }

    for filter_name, filter_kwargs, count in _get_filter_counts():
        if filter_kwargs:
            parameters = dict(json.loads(filter_kwargs))
        else:
            parameters = {}

        parameters["count"] = count
        parameters["were"] = "was" if count == 1 else "were"
        parameters["they"] = "it"  if count == 1 else "they"
        print_err("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    # TODO: Add subsampling in the report template dict now that it's stored in the same table as other filters.
    num_excluded_subsamp = _get_num_excluded_by_subsampling()
    if (args.group_by and args.sequences_per_group) or args.subsample_max_sequences:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print_err(f"\t{num_excluded_subsamp} {'was' if num_excluded_subsamp == 1 else 'were'} dropped because of subsampling criteria{seed_txt}")

    if total_strains_passed == 0:
        empty_results_message = "All samples have been dropped! Check filter rules and metadata file format."
        if args.empty_output_reporting is EmptyOutputReportingMethod.ERROR:
            raise AugurError(empty_results_message)
        elif args.empty_output_reporting is EmptyOutputReportingMethod.WARN:
            print_err(f"WARNING: {empty_results_message}")
        elif args.empty_output_reporting is EmptyOutputReportingMethod.SILENT:
            pass
        else:
            raise ValueError(f"Encountered unhandled --empty-output-reporting method {args.empty_output_reporting!r}")

    print_err(f"{total_strains_passed} {'strain' if total_strains_passed == 1 else 'strains'} passed all filters")


def _get_num_excluded_by_lack_of_metadata():
    """Get number of strains present in other inputs but missing in metadata."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT COUNT(*) AS count
            FROM {constants.INPUT_SOURCE_TABLE}
            WHERE NOT {constants.STRAIN_IN_METADATA_COLUMN}
        """)
        return int(result.fetchone()["count"])


def _get_num_metadata_strains():
    """Returns the number of strains in the original metadata."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT COUNT(*) AS count
            FROM {constants.METADATA_TABLE}
        """)
        return int(result.fetchone()["count"])


def _get_num_excluded_by_subsampling():
    """Returns the number of strains excluded by subsampling."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT COUNT(*) AS count
            FROM {constants.FILTER_REASON_TABLE}
            WHERE {constants.FILTER_REASON_COLUMN} = '{constants.SUBSAMPLE_FILTER_REASON}'
        """)
        return int(result.fetchone()["count"])


# TODO: use _get_valid_strains
def _get_total_strains_passed():
    """Returns the number of strains that pass all filter rules."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT COUNT(*) AS count
            FROM {constants.FILTER_REASON_TABLE}
            WHERE NOT {constants.EXCLUDE_COLUMN} OR {constants.INCLUDE_COLUMN}
        """)
        return int(result.fetchone()["count"])


def _get_filter_counts():
    """
    Returns a tuple for each filter with function name, kwargs, and number of strains included/excluded by it.
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT
                {constants.FILTER_REASON_COLUMN},
                {constants.FILTER_REASON_KWARGS_COLUMN},
                COUNT(*) AS count
            FROM {constants.FILTER_REASON_TABLE}
            WHERE {constants.FILTER_REASON_COLUMN} IS NOT NULL
                AND {constants.FILTER_REASON_COLUMN} != '{constants.SUBSAMPLE_FILTER_REASON}'
            GROUP BY {constants.FILTER_REASON_COLUMN}, {constants.FILTER_REASON_KWARGS_COLUMN}
            ORDER BY {constants.FILTER_REASON_COLUMN}, {constants.FILTER_REASON_KWARGS_COLUMN}
        """)
        for row in result:
            yield (
                str(row[constants.FILTER_REASON_COLUMN]),
                str(row[constants.FILTER_REASON_KWARGS_COLUMN]),
                int(row['count']),
            )
