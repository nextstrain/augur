from collections import defaultdict
import csv
import itertools
import json
import os
import pandas as pd
from tempfile import NamedTemporaryFile

from augur.errors import AugurError
from augur.index import (
    index_sequences,
    index_vcf,
    ID_COLUMN as SEQUENCE_INDEX_ID_COLUMN,
    DELIMITER as SEQUENCE_INDEX_DELIMITER,
)
from augur.io.file import PANDAS_READ_CSV_OPTIONS, open_file
from augur.io.metadata import InvalidDelimiter, Metadata, read_metadata
from augur.io.sequences import read_sequences, write_sequences
from augur.io.print import print_err
from augur.io.vcf import is_vcf as filename_is_vcf, write_vcf
from augur.types import EmptyOutputReportingMethod
from . import include_exclude_rules
from .io import get_useful_metadata_columns, write_metadata_based_outputs
from .include_exclude_rules import apply_filters, construct_filters
from .subsample import subsample


def run(args):
    # Determine whether the sequence index exists or whether should be
    # generated. We need to generate an index if the input sequences are in a
    # VCF, if sequence output has been requested (so we can filter strains by
    # sequences that are present), or if any other sequence-based filters have
    # been requested.
    sequence_strains = None
    sequence_index_path = args.sequence_index
    build_sequence_index = False
    is_vcf = filename_is_vcf(args.sequences)

    # Don't build sequence index with --exclude-all since the only way to add
    # strains back in with this flag are the `--include` or `--include-where`
    # options, so we know we don't need a sequence index to apply any additional
    # filters.
    if sequence_index_path is None and args.sequences and not args.exclude_all:
        build_sequence_index = True

    if build_sequence_index:
        # Generate the sequence index on the fly, for backwards compatibility
        # with older workflows that don't generate the index ahead of time.
        # Create a temporary index using a random filename to avoid collisions
        # between multiple filter commands.
        with NamedTemporaryFile(delete=False) as sequence_index_file:
            sequence_index_path = sequence_index_file.name

        print_err(
            "Note: You did not provide a sequence index, so Augur will generate one.",
            "You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`."
        )

        if is_vcf:
            index_vcf(args.sequences, sequence_index_path)
        else:
            index_sequences(args.sequences, sequence_index_path)

    # Load the sequence index, if a path exists.
    sequence_index = None
    if sequence_index_path:
        sequence_index = pd.read_csv(
            sequence_index_path,
            sep=SEQUENCE_INDEX_DELIMITER,
            index_col=SEQUENCE_INDEX_ID_COLUMN,
            dtype={SEQUENCE_INDEX_ID_COLUMN: "string"},
            **PANDAS_READ_CSV_OPTIONS,
        )

        # Remove temporary index file, if it exists.
        if build_sequence_index:
            os.unlink(sequence_index_path)

        sequence_strains = set(sequence_index.index.values)

    #####################################
    #Filtering steps
    #####################################

    # Load metadata. Metadata are the source of truth for which sequences we
    # want to keep in filtered output.
    valid_strains = set() # TODO: rename this more clearly
    filter_counts = defaultdict(int)

    try:
        metadata_object = Metadata(args.metadata, args.metadata_delimiters, args.metadata_id_columns)
    except InvalidDelimiter:
        raise AugurError(
            f"Could not determine the delimiter of {args.metadata!r}. "
            f"Valid delimiters are: {args.metadata_delimiters!r}. "
            "This can be changed with --metadata-delimiters."
        )
    useful_metadata_columns = get_useful_metadata_columns(args, metadata_object.id_column, metadata_object.columns)

    metadata = read_metadata(
        args.metadata,
        delimiters=[metadata_object.delimiter],
        columns=useful_metadata_columns,
        id_columns=[metadata_object.id_column],
        dtype={col: 'category' for col in useful_metadata_columns},
    )

    duplicate_strains = metadata.index[metadata.index.duplicated()]
    if len(duplicate_strains) > 0:
        raise AugurError(f"The following strains are duplicated in '{args.metadata}':\n" + "\n".join(sorted(duplicate_strains)))

    # FIXME: remove redundant variable from chunking logic
    metadata_strains = set(metadata.index.values)

    # Setup filters.
    exclude_by, include_by = construct_filters(
        args,
        sequence_index,
    )

    # Filter metadata.
    seq_keep, sequences_to_filter, sequences_to_include = apply_filters(
        metadata,
        exclude_by,
        include_by,
    )
    # FIXME: remove redundant variable from chunking logic
    valid_strains = seq_keep

    # Track distinct strains to include, so we can write their
    # corresponding metadata, strains, or sequences later, as needed.
    force_included_strains = {
        record["strain"]
        for record in sequences_to_include
    }

    # Setup logging.
    output_log_context_manager = open_file(args.output_log, "w", newline='')
    output_log_writer = None
    if args.output_log:
        # Log the names of strains that were filtered or force-included, so we
        # can properly account for each strain (e.g., including those that were
        # initially filtered for one reason and then included again for another
        # reason).
        output_log = output_log_context_manager.__enter__()
        output_log_header = ("strain", "filter", "kwargs")
        output_log_writer = csv.DictWriter(
            output_log,
            fieldnames=output_log_header,
            delimiter="\t",
            lineterminator="\n",
        )
        output_log_writer.writeheader()

    # Track reasons for filtered or force-included strains, so we can
    # report total numbers filtered and included at the end. Optionally,
    # write out these reasons to a log file.
    for filtered_strain in itertools.chain(sequences_to_filter, sequences_to_include):
        filter_counts[(filtered_strain["filter"], filtered_strain["kwargs"])] += 1

        # Log the names of strains that were filtered or force-included,
        # so we can properly account for each strain (e.g., including
        # those that were initially filtered for one reason and then
        # included again for another reason).
        if args.output_log:
            output_log_writer.writerow(filtered_strain)

    # Setup grouping. We handle the following major use cases:
    #
    # 1. group by and sequences per group defined -> use the given values by the
    # user to identify the highest priority records from each group in a single
    # pass through the metadata.
    #
    # 2. group by and maximum sequences defined -> use the first pass through
    # the metadata to count the number of records in each group, calculate the
    # sequences per group that satisfies the requested maximum, and use a second
    # pass through the metadata to select that many sequences per group.
    #
    # 3. group by not defined but maximum sequences defined -> use a "dummy"
    # group such that we select at most the requested maximum number of
    # sequences in a single pass through the metadata.
    #
    # Each case relies on a priority queue to track the highest priority records
    # per group. In the best case, we can track these records in a single pass
    # through the metadata. In the worst case, we don't know how many sequences
    # per group to use, so we need to calculate this number after the first pass
    # and use a second pass to add records to the queue.
    group_by = args.group_by or ("_dummy",)

    # Prevent force-included sequences from being included again during
    # subsampling.
    seq_keep = seq_keep - force_included_strains

    if seq_keep and (args.sequences_per_group or args.subsample_max_sequences):
        subsampled_strains = subsample(metadata.loc[list(seq_keep)], args, group_by)
    else:
        subsampled_strains = valid_strains

    num_excluded_subsamp = 0

    # Count and optionally log strains that were not included due to
    # subsampling.
    strains_filtered_by_subsampling = valid_strains - subsampled_strains
    num_excluded_subsamp = len(strains_filtered_by_subsampling)
    if output_log_writer:
        for strain in strains_filtered_by_subsampling:
            output_log_writer.writerow({
                "strain": strain,
                "filter": "subsampling",
                "kwargs": "",
            })

    valid_strains = subsampled_strains

    # Force inclusion of specific strains after filtering and subsampling.
    valid_strains = valid_strains | force_included_strains

    # Write output starting with sequences, if they've been requested. It is
    # possible for the input sequences and sequence index to be out of sync
    # (e.g., the index is a superset of the given sequences input), so we need
    # to update the set of strains to keep based on which strains are actually
    # available.
    if is_vcf:
        if args.output:
            # Get the samples to be deleted, not to keep, for VCF
            dropped_samps = list(sequence_strains - valid_strains)
            write_vcf(args.sequences, args.output, dropped_samps)
    elif args.sequences:
        sequences = read_sequences(args.sequences)

        # If the user requested sequence output, stream to disk all sequences
        # that passed all filters to avoid reading sequences into memory first.
        # Even if we aren't emitting sequences, we track the observed strain
        # names in the sequence file as part of the single pass to allow
        # comparison with the provided sequence index.
        if args.output:
            observed_sequence_strains = set()
            with open_file(args.output, "wt") as output_handle:
                for sequence in sequences:
                    observed_sequence_strains.add(sequence.id)

                    if sequence.id in valid_strains:
                        write_sequences(sequence, output_handle, 'fasta')
        else:
            observed_sequence_strains = {sequence.id for sequence in sequences}

        if sequence_strains != observed_sequence_strains:
            # Warn the user if the expected strains from the sequence index are
            # not a superset of the observed strains.
            if sequence_strains is not None and observed_sequence_strains > sequence_strains:
                print_err(
                    "WARNING: The sequence index is out of sync with the provided sequences.",
                    "Metadata and strain output may not match sequence output."
                )

            # Update the set of available sequence strains.
            sequence_strains = observed_sequence_strains

    if args.output_metadata or args.output_strains:
        write_metadata_based_outputs(args.metadata, args.metadata_delimiters,
                                     args.metadata_id_columns, args.output_metadata,
                                     args.output_strains, valid_strains)

    # Calculate the number of strains that don't exist in either metadata or
    # sequences.
    num_excluded_by_lack_of_metadata = 0
    if sequence_strains:
        num_excluded_by_lack_of_metadata = len(sequence_strains - metadata_strains)


    # Calculate the number of strains passed and filtered.
    total_strains_passed = len(valid_strains)
    total_strains_filtered = len(metadata_strains) + num_excluded_by_lack_of_metadata - total_strains_passed

    print(f"{total_strains_filtered} {'strain was' if total_strains_filtered == 1 else 'strains were'} dropped during filtering")

    if num_excluded_by_lack_of_metadata:
        print(f"\t{num_excluded_by_lack_of_metadata} had no metadata")

    report_template_by_filter_name = {
        include_exclude_rules.filter_by_sequence_index.__name__: "{count} had no sequence data",
        include_exclude_rules.filter_by_exclude_all.__name__: "{count} {were} dropped by `--exclude-all`",
        include_exclude_rules.filter_by_exclude.__name__: "{count} {were} dropped because {they} {were} in {exclude_file}",
        include_exclude_rules.filter_by_exclude_where.__name__: "{count} {were} dropped because of '{exclude_where}'",
        include_exclude_rules.filter_by_query.__name__: "{count} {were} filtered out by the query: \"{query}\"",
        include_exclude_rules.filter_by_ambiguous_date.__name__: "{count} {were} dropped because of their ambiguous date in {ambiguity}",
        include_exclude_rules.filter_by_min_date.__name__: "{count} {were} dropped because {they} {were} earlier than {min_date} or missing a date",
        include_exclude_rules.filter_by_max_date.__name__: "{count} {were} dropped because {they} {were} later than {max_date} or missing a date",
        include_exclude_rules.filter_by_sequence_length.__name__: "{count} {were} dropped because {they} {were} shorter than minimum length of {min_length}bp when only counting standard nucleotide characters A, C, G, or T (case-insensitive)",
        include_exclude_rules.filter_by_non_nucleotide.__name__: "{count} {were} dropped because {they} had non-nucleotide characters",
        include_exclude_rules.skip_group_by_with_ambiguous_year.__name__: "{count} {were} dropped during grouping due to ambiguous year information",
        include_exclude_rules.skip_group_by_with_ambiguous_month.__name__: "{count} {were} dropped during grouping due to ambiguous month information",
        include_exclude_rules.skip_group_by_with_ambiguous_day.__name__: "{count} {were} dropped during grouping due to ambiguous day information",
        include_exclude_rules.force_include_strains.__name__: "{count} {were} added back because {they} {were} in {include_file}",
        include_exclude_rules.force_include_where.__name__: "{count} {were} added back because of '{include_where}'",
    }
    for (filter_name, filter_kwargs), count in filter_counts.items():
        if filter_kwargs:
            parameters = dict(json.loads(filter_kwargs))
        else:
            parameters = {}

        parameters["count"] = count
        parameters["were"] = "was" if count == 1 else "were"
        parameters["they"] = "it"  if count == 1 else "they"
        print("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    if (group_by and args.sequences_per_group) or args.subsample_max_sequences:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print(f"\t{num_excluded_subsamp} {'was' if num_excluded_subsamp == 1 else 'were'} dropped because of subsampling criteria{seed_txt}")

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

    print(f"{total_strains_passed} {'strain' if total_strains_passed == 1 else 'strains'} passed all filters")
