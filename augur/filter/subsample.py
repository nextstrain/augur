from collections import defaultdict
import heapq
import itertools
import uuid
import numpy as np
import pandas as pd
from textwrap import dedent
from typing import Collection, Dict, Iterable, List, Optional, Set, Tuple, Union

from augur.dates import get_year_week, get_year_month_day
from augur.errors import AugurError
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.print import print_err, _n
from . import constants
from .weights_file import WEIGHTS_COLUMN, COLUMN_VALUE_FOR_DEFAULT_WEIGHT, get_default_weight, get_weighted_columns, read_weights_file

Group = Tuple[str, ...]
"""Combination of grouping column values in tuple form."""


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
        A list of metadata (or generated) columns to group records by.

    Returns
    -------
    dict :
        A mapping of strain names to tuples corresponding to the values of the strain's group.

    Examples
    --------
    >>> strains = ["strain1", "strain2"]
    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "2020-01-01", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by = ["region"]
    >>> group_by_strain = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': ('Africa',), 'strain2': ('Europe',)}

    If we group by year or month, these groups are generated from the date
    string.

    >>> group_by = ["year", "month"]
    >>> group_by_strain = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': ('2020', '2020-01'), 'strain2': ('2020', '2020-02')}

    If we omit the grouping columns, the result will group by a dummy column.

    >>> group_by_strain = get_groups_for_subsampling(strains, metadata)
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
    >>> group_by_strain = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': ('2020', '2020-01', 'unknown'), 'strain2': ('2020', '2020-02', 'unknown')}

    We can group metadata without any non-ID columns.

    >>> metadata = pd.DataFrame([{"strain": "strain1"}, {"strain": "strain2"}]).set_index("strain")
    >>> get_groups_for_subsampling(strains, metadata, group_by=('_dummy',))
    {'strain1': ('_dummy',), 'strain2': ('_dummy',)}
    """
    metadata = metadata.loc[list(strains)]
    group_by_strain = {}

    if len(metadata) == 0:
        return group_by_strain

    if not group_by or group_by == ('_dummy',):
        group_by_strain = {strain: ('_dummy',) for strain in strains}
        return group_by_strain

    group_by_set = set(group_by)
    generated_columns_requested = constants.GROUP_BY_GENERATED_COLUMNS & group_by_set

    # If we could not find any requested categories, we cannot complete subsampling.
    if METADATA_DATE_COLUMN not in metadata and group_by_set <= constants.GROUP_BY_GENERATED_COLUMNS:
        raise AugurError(f"The specified group-by categories ({group_by}) were not found. Note that using any of {sorted(constants.GROUP_BY_GENERATED_COLUMNS)} requires a column called {METADATA_DATE_COLUMN!r}.")
    if not group_by_set & (set(metadata.columns) | constants.GROUP_BY_GENERATED_COLUMNS):
        raise AugurError(f"The specified group-by categories ({group_by}) were not found.")

    # Warn/error based on other columns grouped with week.
    if constants.DATE_WEEK_COLUMN in group_by_set:
        if constants.DATE_YEAR_COLUMN in group_by_set:
            print_err(f"WARNING: {constants.DATE_YEAR_COLUMN!r} grouping will be ignored since {constants.DATE_WEEK_COLUMN!r} includes ISO year.")
            group_by.remove(constants.DATE_YEAR_COLUMN)
            group_by_set.remove(constants.DATE_YEAR_COLUMN)
            generated_columns_requested.remove(constants.DATE_YEAR_COLUMN)
        if constants.DATE_MONTH_COLUMN in group_by_set:
            raise AugurError(f"{constants.DATE_MONTH_COLUMN!r} and {constants.DATE_WEEK_COLUMN!r} grouping cannot be used together.")

    if generated_columns_requested:

        if METADATA_DATE_COLUMN not in metadata:
            # Set generated columns to 'unknown'.
            print_err(f"WARNING: A {METADATA_DATE_COLUMN!r} column could not be found to group-by {sorted(generated_columns_requested)}.")
            print_err(f"Filtering by group may behave differently than expected!")
            df_dates = pd.DataFrame({col: 'unknown' for col in constants.GROUP_BY_GENERATED_COLUMNS}, index=metadata.index)
            metadata = pd.concat([metadata, df_dates], axis=1)
        else:
            # Create a DataFrame with year/month/day columns as nullable ints.
            # These columns are prefixed to note temporary usage. They are used
            # to generate other columns, and will be discarded at the end.
            temp_prefix = str(uuid.uuid4())
            temp_date_cols = [f'{temp_prefix}year', f'{temp_prefix}month', f'{temp_prefix}day']
            df_dates = pd.DataFrame(
                metadata[METADATA_DATE_COLUMN].apply(get_year_month_day).tolist(),
                columns=temp_date_cols,
                index=metadata.index,
            )
            for col in temp_date_cols:
                df_dates[col] = pd.to_numeric(df_dates[col], errors='coerce').astype(pd.Int64Dtype())

            # Extend metadata with generated date columns
            # Drop the date column since it should not be used for grouping.
            metadata = pd.concat([metadata.drop(METADATA_DATE_COLUMN, axis=1), df_dates], axis=1)

            # Check again if metadata is empty after dropping ambiguous dates.
            if metadata.empty:
                return group_by_strain

            # Generate columns.
            if constants.DATE_YEAR_COLUMN in generated_columns_requested:
                metadata[constants.DATE_YEAR_COLUMN] = metadata[f'{temp_prefix}year'].astype('string')
            if constants.DATE_MONTH_COLUMN in generated_columns_requested:
                metadata[constants.DATE_MONTH_COLUMN] = (
                    metadata[f'{temp_prefix}year'].astype(str) + '-' +
                    metadata[f'{temp_prefix}month'].astype(str).str.zfill(2)
                )
            if constants.DATE_WEEK_COLUMN in generated_columns_requested:
                # Note that week = (year, week) from the date.isocalendar().
                # Do not combine the raw year with the ISO week number alone,
                # since raw year ≠ ISO year.
                metadata[constants.DATE_WEEK_COLUMN] = metadata.apply(lambda row: get_year_week(
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
        print_err(f"WARNING: Some of the specified group-by categories couldn't be found: {', '.join(unknown_groups)}")
        print_err("Filtering by group may behave differently than expected!")
        for group in unknown_groups:
            metadata[group] = 'unknown'

    # Finally, determine groups.
    group_by_strain = dict(zip(metadata.index, metadata[group_by].apply(tuple, axis=1)))
    return group_by_strain


class PriorityQueue:
    """A priority queue implementation that automatically replaces lower priority
    items in the heap with incoming higher priority items.

    Examples
    --------

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


def get_probabilistic_group_sizes(groups, target_group_size, random_seed=None):
    """Create a dictionary of maximum sizes per group.

    Probabilistically generate varying sizes from a Poisson distribution. Make
    at least the given number of maximum attempts to generate sizes for which
    the total of all sizes is greater than zero.

    Examples
    --------
    Get sizes for two groups with a fractional maximum size. Their total
    size should still be an integer value greater than zero.

    >>> groups = ("2015", "2016")
    >>> seed = 314159
    >>> group_sizes = get_probabilistic_group_sizes(groups, 0.1, random_seed=seed)
    >>> int(sum(group_sizes.values())) > 0
    True

    A subsequent run of this function with the same groups and random seed
    should produce the same group sizes.

    >>> more_group_sizes = get_probabilistic_group_sizes(groups, 0.1, random_seed=seed)
    >>> list(group_sizes.values()) == list(more_group_sizes.values())
    True

    """
    assert target_group_size < 1.0

    # For small fractional maximum sizes, it is possible to randomly select
    # maximum queue sizes that all equal zero. When this happens, filtering
    # fails unexpectedly. We make multiple attempts to create queues with
    # maximum sizes greater than zero for at least one queue.
    random_generator = np.random.default_rng(random_seed)
    total_max_size = 0
    attempts = 0
    max_attempts = 100
    max_sizes_per_group = {}

    while total_max_size == 0 and attempts < max_attempts:
        for group in sorted(groups):
            max_sizes_per_group[group] = random_generator.poisson(target_group_size)

        total_max_size = sum(max_sizes_per_group.values())
        attempts += 1

    return max_sizes_per_group


TARGET_SIZE_COLUMN = '_augur_filter_target_size'
INPUT_SIZE_COLUMN = '_augur_filter_input_size'
OUTPUT_SIZE_COLUMN = '_augur_filter_subsampling_output_size'


def get_weighted_group_sizes(
        records_per_group: Dict[Group, int],
        group_by: List[str],
        weights_file: str,
        target_total_size: int,
        output_sizes_file: Optional[str],
        random_seed: Optional[int],
    ) -> Dict[Group, int]:
    """Return target group sizes based on weights defined in ``weights_file``.
    """
    groups = records_per_group.keys()

    weights = read_weights_file(weights_file)

    weighted_columns = get_weighted_columns(weights_file)

    # Other columns in group_by are considered unweighted.
    unweighted_columns = list(set(group_by) - set(weighted_columns))

    if unweighted_columns:
        # This has the side effect of weighting the values *alongside* (rather
        # than within) each weighted group. After dropping unused groups, adjust
        # weights to ensure equal weighting of unweighted columns *within* each
        # weighted group defined by the weighted columns.
        weights = _add_unweighted_columns(weights, groups, group_by, unweighted_columns)

        weights = _handle_incomplete_weights(weights, weights_file, weighted_columns, group_by, groups)
        weights = _drop_unused_groups(weights, groups, group_by)

        weights = _adjust_weights_for_unweighted_columns(weights, weighted_columns, unweighted_columns)
    else:
        weights = _handle_incomplete_weights(weights, weights_file, weighted_columns, group_by, groups)
        weights = _drop_unused_groups(weights, groups, group_by)

    weights = _calculate_weighted_group_sizes(weights, target_total_size, random_seed)

    # Add columns to summarize the input data
    weights[INPUT_SIZE_COLUMN] = weights.apply(lambda row: records_per_group[tuple(row[group_by].values)], axis=1)
    weights[OUTPUT_SIZE_COLUMN] = weights[[INPUT_SIZE_COLUMN, TARGET_SIZE_COLUMN]].min(axis=1)

    # Warn on any under-sampled groups
    for _, row in weights.iterrows():
        if row[INPUT_SIZE_COLUMN] < row[TARGET_SIZE_COLUMN]:
            sequences = _n('sequence', 'sequences', int(row[TARGET_SIZE_COLUMN]))
            are = _n('is', 'are', int(row[INPUT_SIZE_COLUMN]))
            group = list(f'{col}={value!r}' for col, value in row[group_by].items())
            print_err(f"WARNING: Targeted {row[TARGET_SIZE_COLUMN]} {sequences} for group {group} but only {row[INPUT_SIZE_COLUMN]} {are} available.")

    if output_sizes_file:
        weights.to_csv(output_sizes_file, index=False, sep='\t')

    return dict(zip(weights[group_by].apply(tuple, axis=1), weights[TARGET_SIZE_COLUMN]))


def _add_unweighted_columns(
        weights: pd.DataFrame,
        groups: Iterable[Group],
        group_by: List[str],
        unweighted_columns: List[str],
    ) -> pd.DataFrame:
    """Add the unweighted columns to the weights DataFrame.
    
    This is done by extending the existing weights to the newly created groups.
    """

    # Get unique values for each unweighted column.
    values_for_unweighted_columns = defaultdict(set)
    for group in groups:
        # NOTE: The ordering of entries in `group` corresponds to the column
        # names in `group_by`, but only because `get_groups_for_subsampling`
        # conveniently retains the order. This could be more tightly coupled,
        # but it works.
        column_to_value_map = dict(zip(group_by, group))
        for column in unweighted_columns:
            values_for_unweighted_columns[column].add(column_to_value_map[column])

    # Create a DataFrame for all permutations of values in unweighted columns.
    lists = [sorted(values_for_unweighted_columns[column]) for column in unweighted_columns]
    unweighted_permutations = pd.DataFrame(list(itertools.product(*lists)), columns=unweighted_columns)

    return pd.merge(unweighted_permutations, weights, how='cross')


def _drop_unused_groups(
        weights: pd.DataFrame,
        groups: Collection[Group],
        group_by: List[str],
    ) -> pd.DataFrame:
    """Drop any groups from ``weights`` that don't appear in ``groups``.
    """
    weights.set_index(group_by, inplace=True)

    # Pandas only uses MultiIndex if there is more than one column in the index.
    valid_index: Set[Union[Group, str]]
    if len(group_by) > 1:
        valid_index = set(groups)
    else:
        valid_index = set(group[0] for group in groups)

    extra_groups = set(weights.index) - valid_index
    if extra_groups:
        count = len(extra_groups)
        unit = _n("group", "groups", count)
        print_err(f"NOTE: Skipping {count} {unit} due to lack of entries in metadata.")
        weights = weights[weights.index.isin(valid_index)]

    weights.reset_index(inplace=True)

    return weights


def _adjust_weights_for_unweighted_columns(
        weights: pd.DataFrame,
        weighted_columns: List[str],
        unweighted_columns: Collection[str],
    ) -> pd.DataFrame:
    """Adjust weights for unweighted columns to reflect equal weighting within each weighted group.
    """
    columns = _n('column', 'columns', len(unweighted_columns))
    those = _n('that', 'those', len(unweighted_columns))
    print_err(f"NOTE: Weights were not provided for the {columns} {', '.join(repr(col) for col in unweighted_columns)}. Using equal weights across values in {those} {columns}.")        

    weights_grouped = weights.groupby(weighted_columns)
    weights[WEIGHTS_COLUMN] = weights_grouped[WEIGHTS_COLUMN].transform(lambda x: x / len(x))

    return weights


def _calculate_weighted_group_sizes(
        weights: pd.DataFrame,
        target_total_size: int,
        random_seed: Optional[int],
    ) -> pd.DataFrame:
    """Calculate maximum group sizes based on weights.
    """
    weights[TARGET_SIZE_COLUMN] = pd.Series(weights[WEIGHTS_COLUMN] / weights[WEIGHTS_COLUMN].sum() * target_total_size)

    # Group sizes must be whole numbers. Round probabilistically by adding a
    # random number between [0,1) and truncating the decimal part.
    rng = np.random.default_rng(random_seed)
    weights[TARGET_SIZE_COLUMN] = (weights[TARGET_SIZE_COLUMN].add(pd.Series(rng.random(len(weights))))).astype(int)

    return weights


def _handle_incomplete_weights(
        weights: pd.DataFrame,
        weights_file: str,
        weighted_columns: List[str],
        group_by: List[str],
        groups: Iterable[Group],
    ) -> pd.DataFrame:
    """Handle the case where the weights file does not cover all rows in the metadata.
    """
    missing_groups = set(groups) - set(weights[group_by].apply(tuple, axis=1))

    if not missing_groups:
        return weights

    # Collect the column values that are missing weights.
    missing_values_by_column = defaultdict(set)
    for group in missing_groups:
        # NOTE: The ordering of entries in `group` corresponds to the column
        # names in `group_by`, but only because `get_groups_for_subsampling`
        # conveniently retains the order. This could be more tightly coupled,
        # but it works.
        column_to_value_map = dict(zip(group_by, group))
        for column in weighted_columns:
            missing_values_by_column[column].add(column_to_value_map[column])

    columns_with_values = '\n            - '.join(f'{column!r}: {list(sorted(values))}' for column, values in sorted(missing_values_by_column.items()))

    default_weight = get_default_weight(weights, weighted_columns)

    if not default_weight:
        raise AugurError(dedent(f"""\
            The input metadata contains these values under the following columns that are not covered by {weights_file!r}:
            - {columns_with_values}
            To fix this, either:
            (1) specify weights explicitly - add entries to {weights_file!r} for the values above, or
            (2) specify a default weight - add an entry to {weights_file!r} with the value {COLUMN_VALUE_FOR_DEFAULT_WEIGHT!r} for all columns"""))
    else:
        print_err(dedent(f"""\
            WARNING: The input metadata contains these values under the following columns that are not directly covered by {weights_file!r}:
            - {columns_with_values}
            The default weight of {default_weight} will be used for all groups defined by those values."""))

        missing_weights = pd.DataFrame(sorted(missing_groups), columns=group_by)
        missing_weights[WEIGHTS_COLUMN] = default_weight
        return pd.merge(weights, missing_weights, on=[*group_by, WEIGHTS_COLUMN], how='outer')


def create_queues_by_group(max_sizes_per_group):
    return {group: PriorityQueue(max_size)
            for group, max_size in max_sizes_per_group.items()}


def calculate_sequences_per_group(target_max_value, group_sizes, allow_probabilistic=True):
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
    group_sizes : list of int
        A list with the number of sequences in each requested group.
    allow_probabilistic : bool
        Whether to allow probabilistic subsampling when the number of groups
        exceeds the requested maximum.

    Raises
    ------
    TooManyGroupsError
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
            group_sizes,
        )
    except TooManyGroupsError as error:
        if allow_probabilistic:
            print_err(f"WARNING: {error}")
            sequences_per_group = target_max_value / len(group_sizes)
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
        hypothetical_spg: float, group_sizes: Collection[int],
) -> float:
    # calculate how many sequences we'd keep given a hypothetical spg.
    return sum(
        min(hypothetical_spg, group_count)
        for group_count in group_sizes
    )


def _calculate_sequences_per_group(
        target_max_value: int,
        group_sizes: Collection[int]
) -> int:
    """This is partially inspired by
    https://github.com/python/cpython/blob/3.8/Lib/bisect.py

    This should return the spg such that we don't exceed the requested
    number of samples.

    Parameters
    ----------
    target_max_value : int
        the total number of sequences allowed across all groups
    group_sizes : Collection[int]
        the number of sequences in each group

    Returns
    -------
    int
        maximum number of sequences allowed per group to meet the required maximum total
        sequences allowed

    Examples
    --------
    >>> _calculate_sequences_per_group(4, [4, 2])
    2
    >>> _calculate_sequences_per_group(2, [4, 2])
    1
    >>> _calculate_sequences_per_group(1, [4, 2])
    Traceback (most recent call last):
        ...
    augur.filter.subsample.TooManyGroupsError: Asked to provide at most 1 sequences, but there are 2 groups.
    """

    if len(group_sizes) > target_max_value:
        # we have more groups than sequences we are allowed, which is an
        # error.

        raise TooManyGroupsError(
            "Asked to provide at most {} sequences, but there are {} "
            "groups.".format(target_max_value, len(group_sizes)))

    lo = 1
    hi = target_max_value

    while hi - lo > 2:
        mid = (hi + lo) // 2
        if _calculate_total_sequences(mid, group_sizes) <= target_max_value:
            lo = mid
        else:
            hi = mid

    if _calculate_total_sequences(hi, group_sizes) <= target_max_value:
        return int(hi)
    else:
        return int(lo)
