import json
from augur.errors import AugurError
from augur.io.print import print_err
from augur.types import EmptyOutputReportingMethod
from . import constants, include_exclude_rules


def print_report(args):
    """Print a report of how many strains were dropped and reasoning."""
    # Calculate the number of strains that don't exist in either metadata or
    # sequences.
    num_excluded_by_lack_of_metadata = 0
    if constants.sequence_strains:
        num_excluded_by_lack_of_metadata = len(constants.sequence_strains - constants.metadata_strains)

    # Calculate the number of strains passed and filtered.
    total_strains_passed = len(constants.valid_strains)
    total_strains_filtered = len(constants.metadata_strains) + num_excluded_by_lack_of_metadata - total_strains_passed

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
        include_exclude_rules.force_include_strains.__name__: "{count} {were} added back because {they} {were} in {include_file}",
        include_exclude_rules.force_include_where.__name__: "{count} {were} added back because of '{include_where}'",
    }
    for (filter_name, filter_kwargs), count in constants.filter_counts.items():
        if filter_kwargs:
            parameters = dict(json.loads(filter_kwargs))
        else:
            parameters = {}

        parameters["count"] = count
        parameters["were"] = "was" if count == 1 else "were"
        parameters["they"] = "it"  if count == 1 else "they"
        print_err("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    if (args.group_by and args.sequences_per_group) or args.subsample_max_sequences:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print_err(f"\t{constants.num_excluded_subsamp} {'was' if constants.num_excluded_subsamp == 1 else 'were'} dropped because of subsampling criteria{seed_txt}")

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
