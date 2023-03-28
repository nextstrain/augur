from collections import defaultdict
import csv
import itertools
import numpy as np

from augur.errors import AugurError
from augur.io.metadata import InvalidDelimiter, Metadata, read_metadata
from augur.io.tabular_file import InvalidDelimiter
from . import constants
from .io import cleanup_outputs, get_useful_metadata_columns, read_priority_scores, write_metadata_based_outputs, import_sequence_index, read_and_output_sequences
from .include_exclude_rules import apply_filters, construct_filters
from .report import print_report
from .subsample import PriorityQueue, TooManyGroupsError, calculate_sequences_per_group, create_queues_by_group, get_groups_for_subsampling


def run(args):
    import_sequence_index(args)

    #####################################
    #Filtering steps
    #####################################

    # Setup filters.
    exclude_by, include_by = construct_filters(args)

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
    constants.metadata_strains = set()
    constants.valid_strains = set() # TODO: rename this more clearly
    constants.all_sequences_to_include = set()
    constants.filter_counts = defaultdict(int)

    try:
        metadata_object = Metadata(args.metadata, args.metadata_id_columns, delimiters=args.metadata_delimiters)
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
            set(metadata.index[metadata.index.isin(constants.metadata_strains)])
        )
        if len(duplicate_strains) > 0:
            cleanup_outputs(args)
            raise AugurError(f"The following strains are duplicated in '{args.metadata}':\n" + "\n".join(sorted(duplicate_strains)))

        # Maintain list of all strains seen.
        constants.metadata_strains.update(set(metadata.index.values))

        # Filter metadata.
        seq_keep, sequences_to_filter, sequences_to_include = apply_filters(
            metadata,
            exclude_by,
            include_by,
        )
        constants.valid_strains.update(seq_keep)

        # Track distinct strains to include, so we can write their
        # corresponding metadata, strains, or sequences later, as needed.
        distinct_force_included_strains = {
            record["strain"]
            for record in sequences_to_include
        }
        constants.all_sequences_to_include.update(distinct_force_included_strains)

        # Track reasons for filtered or force-included strains, so we can
        # report total numbers filtered and included at the end. Optionally,
        # write out these reasons to a log file.
        for filtered_strain in itertools.chain(sequences_to_filter, sequences_to_include):
            constants.filter_counts[(filtered_strain["filter"], filtered_strain["kwargs"])] += 1

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
            print(f"Sampling probabilistically at {sequences_per_group:0.4f} sequences per group, meaning it is possible to have more than the requested maximum of {args.subsample_max_sequences} sequences after filtering.")
        else:
            print(f"Sampling at {sequences_per_group} per group.")

        if queues_by_group is None:
            # We know all of the possible groups now from the first pass through
            # the metadata, so we can create queues for all groups at once.
            queues_by_group = create_queues_by_group(
                records_per_group.keys(),
                sequences_per_group,
                random_seed=args.subsample_seed,
            )

        # Make a second pass through the metadata, only considering records that
        # have passed filters.
        metadata_reader = read_metadata(
            args.metadata,
            delimiters=[metadata_object.dialect.delimiter],
            columns=useful_metadata_columns,
            id_columns=[metadata_object.id_column],
            chunk_size=args.metadata_chunk_size,
            dtype="string",
        )
        for metadata in metadata_reader:
            # Recalculate groups for subsampling as we loop through the
            # metadata a second time. TODO: We could store these in memory
            # during the first pass, but we want to minimize overall memory
            # usage at the moment.
            seq_keep = set(metadata.index.values) & constants.valid_strains

            # Prevent force-included strains from being considered in this
            # second pass, as in the first pass.
            seq_keep = seq_keep - constants.all_sequences_to_include

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
    constants.num_excluded_subsamp = 0
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
        strains_filtered_by_subsampling = constants.valid_strains - subsampled_strains
        constants.num_excluded_subsamp = len(strains_filtered_by_subsampling)
        if output_log_writer:
            for strain in strains_filtered_by_subsampling:
                output_log_writer.writerow({
                    "strain": strain,
                    "filter": "subsampling",
                    "kwargs": "",
                })

        constants.valid_strains = subsampled_strains

    read_and_output_sequences(args)

    if args.output_metadata or args.output_strains:
        write_metadata_based_outputs(args.metadata, args.metadata_delimiters,
                                     args.metadata_id_columns, args.output_metadata,
                                     args.output_strains, constants.valid_strains)

    print_report(args)
