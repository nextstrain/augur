from collections import defaultdict
from contextlib import nullcontext
import csv
import itertools
import json
import numpy as np
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
from .io import cleanup_outputs, get_useful_metadata_columns, read_priority_scores, write_metadata_based_outputs
from .include_exclude_rules import apply_filters, construct_filters
from .subsample import PriorityQueue, TooManyGroupsError, calculate_sequences_per_group, get_probabilistic_group_sizes, create_queues_by_group, get_groups_for_subsampling, get_weighted_group_sizes


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

    # Setup filters.
    exclude_by, include_by = construct_filters(
        args,
        sequence_index,
    )

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
    group_by = args.group_by
    sequences_per_group = args.sequences_per_group
    records_per_group = None

    if group_by and args.subsample_max_sequences:
        # In this case, we need two passes through the metadata with the first
        # pass used to count the number of records per group.
        records_per_group = defaultdict(int)
    elif not group_by and args.subsample_max_sequences:
        group_by = ("_dummy",)
        sequences_per_group = args.subsample_max_sequences

    # If we are grouping data, use queues to store the highest priority strains
    # for each group. When no priorities are provided, they will be randomly
    # generated.
    queues_by_group = None
    if group_by:
        # Use user-defined priorities, if possible. Otherwise, setup a
        # corresponding dictionary that returns a random float for each strain.
        if args.priority:
            priorities = read_priority_scores(args.priority)
        else:
            random_generator = np.random.default_rng(args.subsample_seed)
            priorities = defaultdict(random_generator.random)

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

    # Load metadata. Metadata are the source of truth for which sequences we
    # want to keep in filtered output.
    metadata_strains = set()
    valid_strains = set() # TODO: rename this more clearly
    all_sequences_to_include = set()
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

    metadata_reader = read_metadata(
        args.metadata,
        delimiters=[metadata_object.delimiter],
        columns=useful_metadata_columns,
        id_columns=[metadata_object.id_column],
        chunk_size=args.metadata_chunk_size,
        dtype="string",
    )
    for metadata in metadata_reader:
        duplicate_strains = (
            set(metadata.index[metadata.index.duplicated()]) |
            (set(metadata.index) & metadata_strains)
        )
        if len(duplicate_strains) > 0:
            cleanup_outputs(args)
            raise AugurError(f"The following strains are duplicated in '{args.metadata}':\n" + "\n".join(sorted(duplicate_strains)))

        # Maintain list of all strains seen.
        metadata_strains.update(set(metadata.index.values))

        # Filter metadata.
        seq_keep, sequences_to_filter, sequences_to_include = apply_filters(
            metadata,
            exclude_by,
            include_by,
        )
        valid_strains.update(seq_keep)

        # Track distinct strains to include, so we can write their
        # corresponding metadata, strains, or sequences later, as needed.
        distinct_force_included_strains = {
            record["strain"]
            for record in sequences_to_include
        }
        all_sequences_to_include.update(distinct_force_included_strains)

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

        if group_by:
            # Prevent force-included sequences from being included again during
            # subsampling.
            seq_keep = seq_keep - distinct_force_included_strains

            # If grouping, track the highest priority metadata records or
            # count the number of records per group. First, we need to get
            # the groups for the given records.
            group_by_strain = get_groups_for_subsampling(
                seq_keep,
                metadata,
                group_by,
            )

            if args.subsample_max_sequences and records_per_group is not None:
                # Count the number of records per group. We will use this
                # information to calculate the number of sequences per group
                # for the given maximum number of requested sequences.
                for group in group_by_strain.values():
                    records_per_group[group] += 1
            else:
                # Track the highest priority records, when we already
                # know the number of sequences allowed per group.
                if queues_by_group is None:
                    queues_by_group = {}

                for strain in sorted(group_by_strain.keys()):
                    # During this first pass, we do not know all possible
                    # groups will be, so we need to build each group's queue
                    # as we first encounter the group.
                    group = group_by_strain[strain]
                    if group not in queues_by_group:
                        queues_by_group[group] = PriorityQueue(
                            max_size=sequences_per_group,
                        )

                    queues_by_group[group].add(
                        metadata.loc[strain],
                        priorities[strain],
                    )

    # In the worst case, we need to calculate sequences per group from the
    # requested maximum number of sequences and the number of sequences per
    # group. Then, we need to make a second pass through the metadata to find
    # the requested number of records.
    if args.subsample_max_sequences and records_per_group is not None:
        if queues_by_group is None:
            # We know all of the possible groups now from the first pass through
            # the metadata, so we can create queues for all groups at once.
            if args.group_by_weights:
                print_err(f"Sampling with weights defined by {args.group_by_weights}.")
                group_sizes = get_weighted_group_sizes(
                    records_per_group,
                    group_by,
                    args.group_by_weights,
                    args.subsample_max_sequences,
                    args.output_group_by_sizes,
                    args.subsample_seed,
                )
            else:
                # Calculate sequences per group. If there are more groups than maximum
                # sequences requested, sequences per group will be a floating point
                # value and subsampling will be probabilistic.
                try:
                    sequences_per_group, probabilistic_used = calculate_sequences_per_group(
                        args.subsample_max_sequences,
                        records_per_group.values(),
                        args.probabilistic_sampling,
                    )
                except TooManyGroupsError as error:
                    raise AugurError(error)

                if (probabilistic_used):
                    print_err(f"Sampling probabilistically at {sequences_per_group:0.4f} sequences per group, meaning it is possible to have more than the requested maximum of {args.subsample_max_sequences} sequences after filtering.")
                    group_sizes = get_probabilistic_group_sizes(
                        records_per_group.keys(),
                        sequences_per_group,
                        random_seed=args.subsample_seed,
                    )
                else:
                    print_err(f"Sampling at {sequences_per_group} per group.")
                    assert type(sequences_per_group) is int
                    group_sizes = {group: sequences_per_group for group in records_per_group.keys()}
            queues_by_group = create_queues_by_group(group_sizes)

        # Make a second pass through the metadata, only considering records that
        # have passed filters.
        metadata_reader = read_metadata(
            args.metadata,
            delimiters=args.metadata_delimiters,
            columns=useful_metadata_columns,
            id_columns=args.metadata_id_columns,
            chunk_size=args.metadata_chunk_size,
            dtype="string",
        )
        for metadata in metadata_reader:
            # Recalculate groups for subsampling as we loop through the
            # metadata a second time. TODO: We could store these in memory
            # during the first pass, but we want to minimize overall memory
            # usage at the moment.
            seq_keep = set(metadata.index.values) & valid_strains

            # Prevent force-included strains from being considered in this
            # second pass, as in the first pass.
            seq_keep = seq_keep - all_sequences_to_include

            group_by_strain = get_groups_for_subsampling(
                seq_keep,
                metadata,
                group_by,
            )

            for strain in sorted(group_by_strain.keys()):
                group = group_by_strain[strain]
                queues_by_group[group].add(
                    metadata.loc[strain],
                    priorities[strain],
                )

    # If we have any records in queues, we have grouped results and need to
    # stream the highest priority records to the requested outputs.
    num_excluded_subsamp = 0
    if queues_by_group:
        # Populate the set of strains to keep from the records in queues.
        subsampled_strains = set()
        for group, queue in queues_by_group.items():
            records = []
            for record in queue.get_items():
                # Each record is a pandas.Series instance. Track the name of the
                # record, so we can output its sequences later.
                subsampled_strains.add(record.name)

                # Construct a data frame of records to simplify metadata output.
                records.append(record)

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
    valid_strains = valid_strains | all_sequences_to_include

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
        # If the user requested sequence output, stream to disk all sequences
        # that passed all filters to avoid reading sequences into memory first.
        # Even if we aren't emitting sequences, we check for duplicates and
        # track the observed strain names in the sequence file as part of the
        # single pass to allow comparison with the provided sequence index.
        observed_sequence_strains = set()
        duplicates = set()
        with open_file(args.output, "wt") if args.output else nullcontext() as output_handle:
            for sequence in read_sequences(args.sequences):
                if sequence.id in observed_sequence_strains:
                    duplicates.add(sequence.id)

                observed_sequence_strains.add(sequence.id)

                if args.output:
                    if sequence.id in valid_strains:
                        write_sequences(sequence, output_handle, 'fasta')

        if duplicates:
            cleanup_outputs(args)
            raise AugurError(f"The following strains are duplicated in '{args.sequences}':\n" + "\n".join(sorted(duplicates)))

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
    for (filter_name, filter_kwargs), count in filter_counts.items():
        if filter_kwargs:
            parameters = dict(json.loads(filter_kwargs))
        else:
            parameters = {}

        parameters["count"] = count
        parameters["were"] = "was" if count == 1 else "were"
        parameters["they"] = "it"  if count == 1 else "they"
        print_err("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    if (group_by and args.sequences_per_group) or args.subsample_max_sequences:
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
