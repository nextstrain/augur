from collections import defaultdict
import csv
import heapq
import itertools
import json
import numpy as np
import os
import pandas as pd
import sys
import uuid
from tempfile import NamedTemporaryFile
from typing import Collection

from augur.dates import get_iso_year_week
from augur.errors import AugurError
from augur.index import index_sequences, index_vcf
from augur.io.file import open_file
from augur.io.metadata import read_metadata
from augur.io.sequences import read_sequences, write_sequences
from augur.io.vcf import is_vcf as filename_is_vcf, write_vcf
from .io import cleanup_outputs, read_priority_scores
from .include_exclude_rules import apply_filters, construct_filters

from . import GROUP_BY_GENERATED_COLUMNS


SEQUENCE_ONLY_FILTERS = (
    "min_length",
    "non_nucleotide",
)


def get_groups_for_subsampling(strains, metadata, group_by=None):
    """Return a list of groups for each given strain based on the corresponding
    metadata and group by column.

    Parameters
    ----------
    strains : list
        A list of strains to get groups for.
    metadata : pandas.DataFrame
        Metadata to inspect for the given strains.
    group_by : list
        A list of metadata (or calculated) columns to group records by.

    Returns
    -------
    dict :
        A mapping of strain names to tuples corresponding to the values of the strain's group.
    list :
        A list of dictionaries with strains that were skipped from grouping and the reason why (see also: `apply_filters` output).


    >>> strains = ["strain1", "strain2"]
    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "2020-01-01", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by = ["region"]
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': ('Africa',), 'strain2': ('Europe',)}
    >>> skipped_strains
    []

    If we group by year or month, these groups are calculated from the date
    string.

    >>> group_by = ["year", "month"]
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': (2020, (2020, 1)), 'strain2': (2020, (2020, 2))}

    If we omit the grouping columns, the result will group by a dummy column.

    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata)
    >>> group_by_strain
    {'strain1': ('_dummy',), 'strain2': ('_dummy',)}

    If we try to group by columns that don't exist, we get an error.

    >>> group_by = ["missing_column"]
    >>> get_groups_for_subsampling(strains, metadata, group_by)
    Traceback (most recent call last):
      ...
    augur.errors.AugurError: The specified group-by categories (['missing_column']) were not found.

    If we try to group by some columns that exist and some that don't, we allow
    grouping to continue and print a warning message to stderr.

    >>> group_by = ["year", "month", "missing_column"]
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': (2020, (2020, 1), 'unknown'), 'strain2': (2020, (2020, 2), 'unknown')}

    If we group by year month and some records don't have that information in
    their date fields, we should skip those records from the group output and
    track which records were skipped for which reasons.

    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, ["year"])
    >>> group_by_strain
    {'strain2': (2020,)}
    >>> skipped_strains
    [{'strain': 'strain1', 'filter': 'skip_group_by_with_ambiguous_year', 'kwargs': ''}]

    Similarly, if we group by month, we should skip records that don't have
    month information in their date fields.

    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "2020", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, ["month"])
    >>> group_by_strain
    {'strain2': ((2020, 2),)}
    >>> skipped_strains
    [{'strain': 'strain1', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}]

    """
    metadata = metadata.loc[list(strains)]
    group_by_strain = {}
    skipped_strains = []

    if metadata.empty:
        return group_by_strain, skipped_strains

    if not group_by or group_by == ('_dummy',):
        group_by_strain = {strain: ('_dummy',) for strain in strains}
        return group_by_strain, skipped_strains

    group_by_set = set(group_by)
    generated_columns_requested = GROUP_BY_GENERATED_COLUMNS & group_by_set

    # If we could not find any requested categories, we cannot complete subsampling.
    if 'date' not in metadata and group_by_set <= GROUP_BY_GENERATED_COLUMNS:
        raise AugurError(f"The specified group-by categories ({group_by}) were not found. Note that using any of {sorted(GROUP_BY_GENERATED_COLUMNS)} requires a column called 'date'.")
    if not group_by_set & (set(metadata.columns) | GROUP_BY_GENERATED_COLUMNS):
        raise AugurError(f"The specified group-by categories ({group_by}) were not found.")

    # Warn/error based on other columns grouped with 'week'.
    if 'week' in group_by_set:
        if 'year' in group_by_set:
            print(f"WARNING: 'year' grouping will be ignored since 'week' includes ISO year.", file=sys.stderr)
            group_by.remove('year')
            group_by_set.remove('year')
            generated_columns_requested.remove('year')
        if 'month' in group_by_set:
            raise AugurError("'month' and 'week' grouping cannot be used together.")

    if generated_columns_requested:

        for col in sorted(generated_columns_requested):
            if col in metadata.columns:
                print(f"WARNING: `--group-by {col}` uses a generated {col} value from the 'date' column. The custom '{col}' column in the metadata is ignored for grouping purposes.", file=sys.stderr)
                metadata.drop(col, axis=1, inplace=True)

        if 'date' not in metadata:
            # Set generated columns to 'unknown'.
            print(f"WARNING: A 'date' column could not be found to group-by {sorted(generated_columns_requested)}.", file=sys.stderr)
            print(f"Filtering by group may behave differently than expected!", file=sys.stderr)
            df_dates = pd.DataFrame({col: 'unknown' for col in GROUP_BY_GENERATED_COLUMNS}, index=metadata.index)
            metadata = pd.concat([metadata, df_dates], axis=1)
        else:
            # Create a DataFrame with year/month/day columns as nullable ints.
            # These columns are prefixed to note temporary usage. They are used
            # to generate other columns, and will be discarded at the end.
            temp_prefix = str(uuid.uuid4())
            temp_date_cols = [f'{temp_prefix}year', f'{temp_prefix}month', f'{temp_prefix}day']
            df_dates = metadata['date'].str.split('-', n=2, expand=True)
            df_dates = df_dates.set_axis(temp_date_cols[:len(df_dates.columns)], axis=1)
            missing_date_cols = set(temp_date_cols) - set(df_dates.columns)
            for col in missing_date_cols:
                df_dates[col] = pd.NA
            for col in temp_date_cols:
                df_dates[col] = pd.to_numeric(df_dates[col], errors='coerce').astype(pd.Int64Dtype())

            # Extend metadata with generated date columns
            # Drop the 'date' column since it should not be used for grouping.
            metadata = pd.concat([metadata.drop('date', axis=1), df_dates], axis=1)

            ambiguous_date_strains = list(_get_ambiguous_date_skipped_strains(
                metadata,
                temp_prefix,
                generated_columns_requested
            ))
            metadata.drop([record['strain'] for record in ambiguous_date_strains], inplace=True)
            skipped_strains.extend(ambiguous_date_strains)

            # Check again if metadata is empty after dropping ambiguous dates.
            if metadata.empty:
                return group_by_strain, skipped_strains

            # Generate columns.
            if 'year' in generated_columns_requested:
                metadata['year'] = metadata[f'{temp_prefix}year']
            if 'month' in generated_columns_requested:
                metadata['month'] = list(zip(
                    metadata[f'{temp_prefix}year'],
                    metadata[f'{temp_prefix}month']
                ))
            if 'week' in generated_columns_requested:
                # Note that week = (year, week) from the date.isocalendar().
                # Do not combine the raw year with the ISO week number alone,
                # since raw year ≠ ISO year.
                metadata['week'] = metadata.apply(lambda row: get_iso_year_week(
                    row[f'{temp_prefix}year'],
                    row[f'{temp_prefix}month'],
                    row[f'{temp_prefix}day']
                    ), axis=1
                )

            # Drop the internally used columns.
            for col in temp_date_cols:
                metadata.drop(col, axis=1, inplace=True)

    unknown_groups = group_by_set - set(metadata.columns)
    if unknown_groups:
        print(f"WARNING: Some of the specified group-by categories couldn't be found: {', '.join(unknown_groups)}", file=sys.stderr)
        print("Filtering by group may behave differently than expected!", file=sys.stderr)
        for group in unknown_groups:
            metadata[group] = 'unknown'

    # Finally, determine groups.
    group_by_strain = dict(zip(metadata.index, metadata[group_by].apply(tuple, axis=1)))
    return group_by_strain, skipped_strains


def _get_ambiguous_date_skipped_strains(
        metadata, temp_prefix, generated_columns_requested):
    """Get strains skipped due to date ambiguity.

    Each value is a dictionary with keys:
    - `strain`: strain name
    - `filter`: filter reason. Used for the final report output.
    - `kwargs`: Empty string since filter reason does not represent a function.
    """
    # Don't yield the same strain twice.
    already_skipped_strains = set()

    if generated_columns_requested:
        # Skip ambiguous years.
        df_skip = metadata[metadata[f'{temp_prefix}year'].isnull()]
        for strain in df_skip.index:
            if strain not in already_skipped_strains:
                yield {
                    "strain": strain,
                    "filter": "skip_group_by_with_ambiguous_year",
                    "kwargs": "",
                }
        already_skipped_strains.update(df_skip.index)
    if 'month' in generated_columns_requested or 'week' in generated_columns_requested:
        # Skip ambiguous months.
        df_skip = metadata[metadata[f'{temp_prefix}month'].isnull()]
        for strain in df_skip.index:
            if strain not in already_skipped_strains:
                yield {
                    "strain": strain,
                    "filter": "skip_group_by_with_ambiguous_month",
                    "kwargs": "",
                }
        already_skipped_strains.update(df_skip.index)
    if 'week' in generated_columns_requested:
        # Skip ambiguous days.
        df_skip = metadata[metadata[f'{temp_prefix}day'].isnull()]
        for strain in df_skip.index:
            if strain not in already_skipped_strains:
                yield {
                    "strain": strain,
                    "filter": "skip_group_by_with_ambiguous_day",
                    "kwargs": "",
                }
        # TODO: uncomment if another filter reason is ever added.
        # already_skipped_strains.update(df_skip.index)


class PriorityQueue:
    """A priority queue implementation that automatically replaces lower priority
    items in the heap with incoming higher priority items.

    Add a single record to a heap with a maximum of 2 records.

    >>> queue = PriorityQueue(max_size=2)
    >>> queue.add({"strain": "strain1"}, 0.5)
    1

    Add another record with a higher priority. The queue should be at its maximum
    size.

    >>> queue.add({"strain": "strain2"}, 1.0)
    2
    >>> queue.heap
    [(0.5, 0, {'strain': 'strain1'}), (1.0, 1, {'strain': 'strain2'})]
    >>> list(queue.get_items())
    [{'strain': 'strain1'}, {'strain': 'strain2'}]

    Add a higher priority record that causes the queue to exceed its maximum
    size. The resulting queue should contain the two highest priority records
    after the lowest priority record is removed.

    >>> queue.add({"strain": "strain3"}, 2.0)
    2
    >>> list(queue.get_items())
    [{'strain': 'strain2'}, {'strain': 'strain3'}]

    Add a record with the same priority as another record, forcing the duplicate
    to be resolved by removing the oldest entry.

    >>> queue.add({"strain": "strain4"}, 1.0)
    2
    >>> list(queue.get_items())
    [{'strain': 'strain4'}, {'strain': 'strain3'}]

    """
    def __init__(self, max_size):
        """Create a fixed size heap (priority queue)

        """
        self.max_size = max_size
        self.heap = []
        self.counter = itertools.count()

    def add(self, item, priority):
        """Add an item to the queue with a given priority.

        If adding the item causes the queue to exceed its maximum size, replace
        the lowest priority item with the given item. The queue stores items
        with an additional heap id value (a count) to resolve ties between items
        with equal priority (favoring the most recently added item).

        """
        heap_id = next(self.counter)

        if len(self.heap) >= self.max_size:
            heapq.heappushpop(self.heap, (priority, heap_id, item))
        else:
            heapq.heappush(self.heap, (priority, heap_id, item))

        return len(self.heap)

    def get_items(self):
        """Return each item in the queue in order.

        Yields
        ------
        Any
            Item stored in the queue.

        """
        for priority, heap_id, item in self.heap:
            yield item


def create_queues_by_group(groups, max_size, max_attempts=100, random_seed=None):
    """Create a dictionary of priority queues per group for the given maximum size.

    When the maximum size is fractional, probabilistically sample the maximum
    size from a Poisson distribution. Make at least the given number of maximum
    attempts to create queues for which the sum of their maximum sizes is
    greater than zero.

    Create queues for two groups with a fixed maximum size.

    >>> groups = ("2015", "2016")
    >>> queues = create_queues_by_group(groups, 2)
    >>> sum(queue.max_size for queue in queues.values())
    4

    Create queues for two groups with a fractional maximum size. Their total max
    size should still be an integer value greater than zero.

    >>> seed = 314159
    >>> queues = create_queues_by_group(groups, 0.1, random_seed=seed)
    >>> int(sum(queue.max_size for queue in queues.values())) > 0
    True

    A subsequent run of this function with the same groups and random seed
    should produce the same queues and queue sizes.

    >>> more_queues = create_queues_by_group(groups, 0.1, random_seed=seed)
    >>> [queue.max_size for queue in queues.values()] == [queue.max_size for queue in more_queues.values()]
    True

    """
    queues_by_group = {}
    total_max_size = 0
    attempts = 0

    if max_size < 1.0:
        random_generator = np.random.default_rng(random_seed)

    # For small fractional maximum sizes, it is possible to randomly select
    # maximum queue sizes that all equal zero. When this happens, filtering
    # fails unexpectedly. We make multiple attempts to create queues with
    # maximum sizes greater than zero for at least one queue.
    while total_max_size == 0 and attempts < max_attempts:
        for group in sorted(groups):
            if max_size < 1.0:
                queue_max_size = random_generator.poisson(max_size)
            else:
                queue_max_size = max_size

            queues_by_group[group] = PriorityQueue(queue_max_size)

        total_max_size = sum(queue.max_size for queue in queues_by_group.values())
        attempts += 1

    return queues_by_group


def validate_arguments(args):
    """Validate arguments and return a boolean representing whether all validation
    rules succeeded.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from argparse

    Returns
    -------
    bool :
        Validation succeeded.

    """
    # Don't allow sequence output when no sequence input is provided.
    if args.output and not args.sequences:
        print(
            "ERROR: You need to provide sequences to output sequences.",
            file=sys.stderr)
        return False

    # Confirm that at least one output was requested.
    if not any((args.output, args.output_metadata, args.output_strains)):
        print(
            "ERROR: You need to select at least one output.",
            file=sys.stderr)
        return False

    # Don't allow filtering on sequence-based information, if no sequences or
    # sequence index is provided.
    if not args.sequences and not args.sequence_index and any(getattr(args, arg) for arg in SEQUENCE_ONLY_FILTERS):
        print(
            "ERROR: You need to provide a sequence index or sequences to filter on sequence-specific information.",
            file=sys.stderr)
        return False

    # Set flags if VCF
    is_vcf = filename_is_vcf(args.sequences)

    # Confirm that vcftools is installed.
    if is_vcf:
        from shutil import which
        if which("vcftools") is None:
            print("ERROR: 'vcftools' is not installed! This is required for VCF data. "
                  "Please see the augur install instructions to install it.",
                  file=sys.stderr)
            return False

    # If user requested grouping, confirm that other required inputs are provided, too.
    if args.group_by and not any((args.sequences_per_group, args.subsample_max_sequences)):
        print(
            "ERROR: You must specify a number of sequences per group or maximum sequences to subsample.",
            file=sys.stderr
        )
        return False

    return True


def run(args):
    # Validate arguments before attempting any I/O.
    if not validate_arguments(args):
        return 1

    # Determine whether the sequence index exists or whether should be
    # generated. We need to generate an index if the input sequences are in a
    # VCF, if sequence output has been requested (so we can filter strains by
    # sequences that are present), or if any other sequence-based filters have
    # been requested.
    sequence_strains = None
    sequence_index_path = args.sequence_index
    build_sequence_index = False
    is_vcf = filename_is_vcf(args.sequences)

    if sequence_index_path is None and args.sequences and not args.exclude_all:
        build_sequence_index = True

    if build_sequence_index:
        # Generate the sequence index on the fly, for backwards compatibility
        # with older workflows that don't generate the index ahead of time.
        # Create a temporary index using a random filename to avoid collisions
        # between multiple filter commands.
        with NamedTemporaryFile(delete=False) as sequence_index_file:
            sequence_index_path = sequence_index_file.name

        print(
            "Note: You did not provide a sequence index, so Augur will generate one.",
            "You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.",
            file=sys.stderr
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
            sep="\t",
            index_col="strain",
        )

        # Remove temporary index file, if it exists.
        if build_sequence_index:
            os.unlink(sequence_index_path)

        # Calculate summary statistics needed for filtering.
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

    # Setup metadata output. We track whether any records have been written to
    # disk yet through the following variables, to control whether we write the
    # metadata's header and open a new file for writing.
    metadata_header = True
    metadata_mode = "w"

    # Setup strain output.
    if args.output_strains:
        output_strains = open(args.output_strains, "w")

    # Setup logging.
    output_log_writer = None
    if args.output_log:
        # Log the names of strains that were filtered or force-included, so we
        # can properly account for each strain (e.g., including those that were
        # initially filtered for one reason and then included again for another
        # reason).
        output_log = open(args.output_log, "w", newline='')
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

    metadata_reader = read_metadata(
        args.metadata,
        id_columns=args.metadata_id_columns,
        chunk_size=args.metadata_chunk_size,
    )
    for metadata in metadata_reader:
        duplicate_strains = (
            set(metadata.index[metadata.index.duplicated()]) |
            set(metadata.index[metadata.index.isin(metadata_strains)])
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
            group_by_strain, skipped_strains = get_groups_for_subsampling(
                seq_keep,
                metadata,
                group_by,
            )

            # Track strains skipped during grouping, so users know why those
            # strains were excluded from the analysis.
            for skipped_strain in skipped_strains:
                filter_counts[(skipped_strain["filter"], skipped_strain["kwargs"])] += 1
                valid_strains.remove(skipped_strain["strain"])

                if args.output_log:
                    output_log_writer.writerow(skipped_strain)

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

        # Always write out strains that are force-included. Additionally, if
        # we are not grouping, write out metadata and strains that passed
        # filters so far.
        force_included_strains_to_write = distinct_force_included_strains
        if not group_by:
            force_included_strains_to_write = force_included_strains_to_write | seq_keep

        if args.output_metadata:
            # TODO: wrap logic to write metadata into its own function
            metadata.loc[list(force_included_strains_to_write)].to_csv(
                args.output_metadata,
                sep="\t",
                header=metadata_header,
                mode=metadata_mode,
            )
            metadata_header = False
            metadata_mode = "a"

        if args.output_strains:
            # TODO: Output strains will no longer be ordered. This is a
            # small breaking change.
            for strain in force_included_strains_to_write:
                output_strains.write(f"{strain}\n")

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
            print(f"ERROR: {error}", file=sys.stderr)
            sys.exit(1)

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
            id_columns=args.metadata_id_columns,
            chunk_size=args.metadata_chunk_size,
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

            group_by_strain, skipped_strains = get_groups_for_subsampling(
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

                if args.output_strains:
                    # TODO: Output strains will no longer be ordered. This is a
                    # small breaking change.
                    output_strains.write(f"{record.name}\n")

            # Write records to metadata output, if requested.
            if args.output_metadata and len(records) > 0:
                records = pd.DataFrame(records)
                records.to_csv(
                    args.output_metadata,
                    sep="\t",
                    header=metadata_header,
                    mode=metadata_mode,
                )
                metadata_header = False
                metadata_mode = "a"

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
                print(
                    "WARNING: The sequence index is out of sync with the provided sequences.",
                    "Metadata and strain output may not match sequence output.",
                    file=sys.stderr
                )

            # Update the set of available sequence strains.
            sequence_strains = observed_sequence_strains

    # Calculate the number of strains that don't exist in either metadata or
    # sequences.
    num_excluded_by_lack_of_metadata = 0
    if sequence_strains:
        # Update strains to keep based on available sequence data. This prevents
        # writing out strain lists or metadata for strains that have no
        # sequences.
        valid_strains = valid_strains & sequence_strains

        num_excluded_by_lack_of_metadata = len(sequence_strains - metadata_strains)

    if args.output_strains:
        output_strains.close()

    # Calculate the number of strains passed and filtered.
    total_strains_passed = len(valid_strains)
    total_strains_filtered = len(metadata_strains) + num_excluded_by_lack_of_metadata - total_strains_passed

    print(f"{total_strains_filtered} strains were dropped during filtering")

    if num_excluded_by_lack_of_metadata:
        print(f"\t{num_excluded_by_lack_of_metadata} had no metadata")

    report_template_by_filter_name = {
        "filter_by_sequence_index": "{count} had no sequence data",
        "filter_by_exclude_all": "{count} of these were dropped by `--exclude-all`",
        "filter_by_exclude": "{count} of these were dropped because they were in {exclude_file}",
        "filter_by_exclude_where": "{count} of these were dropped because of '{exclude_where}'",
        "filter_by_query": "{count} of these were filtered out by the query: \"{query}\"",
        "filter_by_ambiguous_date": "{count} of these were dropped because of their ambiguous date in {ambiguity}",
        "filter_by_min_date": "{count} of these were dropped because they were earlier than {min_date} or missing a date",
        "filter_by_max_date": "{count} of these were dropped because they were later than {max_date} or missing a date",
        "filter_by_sequence_length": "{count} of these were dropped because they were shorter than minimum length of {min_length}bp",
        "filter_by_non_nucleotide": "{count} of these were dropped because they had non-nucleotide characters",
        "skip_group_by_with_ambiguous_year": "{count} were dropped during grouping due to ambiguous year information",
        "skip_group_by_with_ambiguous_month": "{count} were dropped during grouping due to ambiguous month information",
        "skip_group_by_with_ambiguous_day": "{count} were dropped during grouping due to ambiguous day information",
        "force_include_strains": "{count} strains were added back because they were in {include_file}",
        "force_include_where": "{count} sequences were added back because of '{include_where}'",
    }
    for (filter_name, filter_kwargs), count in filter_counts.items():
        if filter_kwargs:
            parameters = dict(json.loads(filter_kwargs))
        else:
            parameters = {}

        parameters["count"] = count
        print("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    if (group_by and args.sequences_per_group) or args.subsample_max_sequences:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print("\t%i of these were dropped because of subsampling criteria%s" % (num_excluded_subsamp, seed_txt))

    if total_strains_passed == 0:
        print("ERROR: All samples have been dropped! Check filter rules and metadata file format.", file=sys.stderr)
        return 1

    print(f"{total_strains_passed} strains passed all filters")


def calculate_sequences_per_group(target_max_value, counts_per_group, allow_probabilistic=True):
    """Calculate the number of sequences per group for a given maximum number of
    sequences to be returned and the number of sequences in each requested
    group. Optionally, allow the result to be probabilistic such that the mean
    result of a Poisson process achieves the calculated sequences per group for
    the given maximum.

    Parameters
    ----------
    target_max_value : int
        Maximum number of sequences to return by subsampling at some calculated
        number of sequences per group for the given counts per group.
    counts_per_group : list[int]
        A list with the number of sequences in each requested group.
    allow_probabilistic : bool
        Whether to allow probabilistic subsampling when the number of groups
        exceeds the requested maximum.

    Raises
    ------
    TooManyGroupsError :
        When there are more groups than sequences per group and probabilistic
        subsampling is not allowed.

    Returns
    -------
    int or float :
        Number of sequences per group.
    bool :
        Whether probabilistic subsampling was used.

    """
    probabilistic_used = False

    try:
        sequences_per_group = _calculate_sequences_per_group(
            target_max_value,
            counts_per_group,
        )
    except TooManyGroupsError as error:
        if allow_probabilistic:
            print(f"WARNING: {error}", file=sys.stderr)
            sequences_per_group = _calculate_fractional_sequences_per_group(
                target_max_value,
                counts_per_group,
            )
            probabilistic_used = True
        else:
            raise error

    return sequences_per_group, probabilistic_used


class TooManyGroupsError(ValueError):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)


def _calculate_total_sequences(
        hypothetical_spg: float, sequence_lengths: Collection[int],
) -> float:
    # calculate how many sequences we'd keep given a hypothetical spg.
    return sum(
        min(hypothetical_spg, sequence_length)
        for sequence_length in sequence_lengths
    )


def _calculate_sequences_per_group(
        target_max_value: int,
        sequence_lengths: Collection[int]
) -> int:
    """This is partially inspired by
    https://github.com/python/cpython/blob/3.8/Lib/bisect.py

    This should return the spg such that we don't exceed the requested
    number of samples.

    Parameters
    ----------
    target_max_value : int
        the total number of sequences allowed across all groups
    sequence_lengths : Collection[int]
        the number of sequences in each group

    Returns
    -------
    int
        maximum number of sequences allowed per group to meet the required maximum total
        sequences allowed

    >>> _calculate_sequences_per_group(4, [4, 2])
    2
    >>> _calculate_sequences_per_group(2, [4, 2])
    1
    >>> _calculate_sequences_per_group(1, [4, 2])
    Traceback (most recent call last):
        ...
    augur.filter.TooManyGroupsError: Asked to provide at most 1 sequences, but there are 2 groups.
    """

    if len(sequence_lengths) > target_max_value:
        # we have more groups than sequences we are allowed, which is an
        # error.

        raise TooManyGroupsError(
            "Asked to provide at most {} sequences, but there are {} "
            "groups.".format(target_max_value, len(sequence_lengths)))

    lo = 1
    hi = target_max_value

    while hi - lo > 2:
        mid = (hi + lo) // 2
        if _calculate_total_sequences(mid, sequence_lengths) <= target_max_value:
            lo = mid
        else:
            hi = mid

    if _calculate_total_sequences(hi, sequence_lengths) <= target_max_value:
        return int(hi)
    else:
        return int(lo)


def _calculate_fractional_sequences_per_group(
        target_max_value: int,
        sequence_lengths: Collection[int]
) -> float:
    """Returns the fractional sequences per group for the given list of group
    sequences such that the total doesn't exceed the requested number of
    samples.

    Parameters
    ----------
    target_max_value : int
        the total number of sequences allowed across all groups
    sequence_lengths : Collection[int]
        the number of sequences in each group

    Returns
    -------
    float
        fractional maximum number of sequences allowed per group to meet the
        required maximum total sequences allowed

    >>> np.around(_calculate_fractional_sequences_per_group(4, [4, 2]), 4)
    1.9375
    >>> np.around(_calculate_fractional_sequences_per_group(2, [4, 2]), 4)
    0.9688

    Unlike the integer-based version of this function, the fractional version
    can accept a maximum number of sequences that exceeds the number of groups.
    In this case, the function returns a fraction that can be used downstream,
    for example with Poisson sampling.

    >>> np.around(_calculate_fractional_sequences_per_group(1, [4, 2]), 4)
    0.4844
    """
    lo = 1e-5
    hi = target_max_value

    while (hi / lo) > 1.1:
        mid = (lo + hi) / 2
        if _calculate_total_sequences(mid, sequence_lengths) <= target_max_value:
            lo = mid
        else:
            hi = mid

    return (lo + hi) / 2