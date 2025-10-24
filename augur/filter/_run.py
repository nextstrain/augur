from collections import defaultdict
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
from augur.io.sequences import read_sequence_ids, subset_fasta, is_vcf as filename_is_vcf, write_vcf
from augur.io.print import print_debug, print_err, _n
from augur.types import EmptyOutputReportingMethod
from . import include_exclude_rules
from .io import cleanup_outputs, get_useful_metadata_columns, read_priority_scores, write_output_metadata
from .include_exclude_rules import apply_filters, construct_filters
from .subsample import PriorityQueue, TooManyGroupsError, calculate_sequences_per_group, get_probabilistic_group_sizes, create_queues_by_group, get_groups_for_subsampling, get_weighted_group_sizes
from .validate_arguments import SEQUENCE_ONLY_FILTERS


def run(args):
    # Determine whether the sequence index exists or whether should be
    # generated.
    sequence_strains = None
    sequence_index_path = args.sequence_index
    sequence_index_ids = None
    build_sequence_index = False
    is_vcf = filename_is_vcf(args.sequences)

    # Build the sequence index if an input sequence file is provided and
    # sequence-based filters will be used. It can be skipped with --exclude-all
    # since the only way to add strains back in with this flag are the
    # `--include` or `--include-where` options.
    if (sequence_index_path is None and
        args.sequences and
        any(getattr(args, arg) for arg in SEQUENCE_ONLY_FILTERS) and
        not args.exclude_all):
        build_sequence_index = True

    if build_sequence_index:
        print_debug(f"Building sequence index for {args.sequences!r}…")
        # Generate the sequence index on the fly for workflows that don't do
        # this separately. Create a temporary index using a random filename to
        # avoid collisions between multiple filter commands.
        with NamedTemporaryFile(delete=False) as sequence_index_file:
            sequence_index_path = sequence_index_file.name

        if is_vcf:
            index_vcf(args.sequences, sequence_index_path)
        else:
            index_sequences(args.sequences, sequence_index_path)

    # Load the sequence index, if a path exists.
    sequence_index = None
    if sequence_index_path:
        if args.sequence_index:
            print_debug(f"Reading sequence index from {sequence_index_path!r}…")
        else:
            print_debug("Reading sequence index from temporary file…")
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

        sequence_index_ids = set(sequence_index.index.values)

    # Check ids for duplicates and compare against sequence index.
    if args.skip_checks:
        print_debug(f"Skipping first pass of sequence file due to --skip-checks.")
    elif args.sequences:
        print_debug(f"Reading sequences from {args.sequences!r}…")

        try:
            sequence_strains = read_sequence_ids(args.sequences, args.nthreads)
        except AugurError as e:
            cleanup_outputs(args)
            raise e

        # Warn the user if the expected strains from the sequence index are
        # not a superset of the observed strains.
        if args.sequence_index and sequence_strains > sequence_index_ids:
            print_err(
                "WARNING: The sequence index is out of sync with the provided sequences. ",
                "Sequence-based filters may drop sequences that are present in --sequences that would have otherwise passed the filters."
            )

    # Use ids from sequence index if no sequence file is provided.
    if sequence_strains is None and sequence_index_ids:
        sequence_strains = sequence_index_ids

    #####################################
    #Filtering steps
    #####################################

    # Setup filters.
    exclude_by, include_by = construct_filters(
        args,
        sequence_index,
        sequence_strains,
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

    print_debug(f"Reading metadata from {args.metadata!r}…")
    metadata_reader = read_metadata(
        args.metadata,
        delimiters=[metadata_object.delimiter],
        columns=useful_metadata_columns,
        id_columns=[metadata_object.id_column],
        keep_id_as_column=True,
        chunk_size=args.metadata_chunk_size,
        dtype="string",
    )
    for metadata in metadata_reader:
        if len(metadata.loc[metadata.index == '']):
            cleanup_outputs(args)
            raise AugurError(f"Found rows with empty values in id column {metadata.index.name!r} in {args.metadata!r}\n" + \
                             "Please remove the rows with empty ids or use a different id column via --metadata-id-columns.")

        duplicate_strains = (
            set(metadata.index[metadata.index.duplicated()]) |
            (set(metadata.index) & metadata_strains)
        )
        if len(duplicate_strains) > 0:
            cleanup_outputs(args)
            raise AugurError(f"The following strains are duplicated in '{args.metadata}':\n" + "\n".join(repr(x) for x in sorted(duplicate_strains)))

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
                        strain,
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

        print_debug(f"Reading metadata from {args.metadata!r}…")
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
                    strain,
                    priorities[strain],
                )

    # If we have any records in queues, we have grouped results and need to
    # stream the highest priority records to the requested outputs.
    num_excluded_subsamp = 0
    if queues_by_group:
        # Populate the set of strains to keep from the records in queues.
        subsampled_strains = set()
        for group, queue in queues_by_group.items():
            for strain in queue.get_items():
                subsampled_strains.add(strain)

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

    # Write requested outputs.

    # Write a strains file if explicitly requested or if sequence output is requested.
    strains_file = None
    if args.output_strains:
        strains_file = args.output_strains
    elif args.output_sequences:
        strains_file = NamedTemporaryFile(delete=False).name

    if strains_file is not None:
        print_debug(f"Writing strains to {strains_file!r}…")
        with open(strains_file, "w") as f:
            for strain in valid_strains:
                f.write(f"{strain}\n")

    if args.output_sequences:
        print_debug(f"Reading sequences from {args.sequences!r} and writing to {args.output_sequences!r}…")
        if is_vcf:
            # Get the samples to be deleted, not to keep, for VCF
            dropped_samps = list(sequence_strains - valid_strains)
            write_vcf(args.sequences, args.output_sequences, dropped_samps)
        else:
            subset_fasta(args.sequences, args.output_sequences, strains_file, args.nthreads)
            if not args.output_strains:
                os.remove(strains_file)

    if args.output_metadata:
        print_debug(f"Reading metadata from {args.metadata!r} and writing to {args.output_metadata!r}…")
        write_output_metadata(args.metadata, args.metadata_delimiters,
                              args.metadata_id_columns, args.output_metadata,
                              valid_strains)

    # Calculate the number of strains that don't exist in either metadata or
    # sequences.
    num_excluded_by_lack_of_metadata = 0
    if sequence_strains:
        num_excluded_by_lack_of_metadata = len(sequence_strains - metadata_strains)


    # Calculate the number of strains passed and filtered.
    total_strains_passed = len(valid_strains)
    total_strains_filtered = len(metadata_strains) + num_excluded_by_lack_of_metadata - total_strains_passed

    print_err(f"{total_strains_filtered} {_n('strain was', 'strains were', total_strains_filtered)} dropped during filtering")

    if num_excluded_by_lack_of_metadata:
        print_err(f"\t{num_excluded_by_lack_of_metadata} had no metadata")

    report_template_by_filter_name = {
        include_exclude_rules.filter_by_sequence_ids.__name__: "{count} had no sequence data",
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
        parameters["were"] = _n("was", "were", count)
        parameters["they"] = _n("it", "they", count)
        print_err("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    if (group_by and args.sequences_per_group) or args.subsample_max_sequences:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print_err(f"\t{num_excluded_subsamp} {_n('was', 'were', num_excluded_subsamp)} dropped because of subsampling criteria{seed_txt}")

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

    print_err(f"{total_strains_passed} {_n('strain', 'strains', total_strains_passed)} passed all filters")
