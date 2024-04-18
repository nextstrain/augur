import uuid
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Collection

from augur.dates import get_iso_year_week
from augur.errors import AugurError
from augur.filter.io import read_priority_scores
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.print import print_err
from . import constants


def subsample(metadata, args, group_by):

    # Use user-defined priorities, if possible. Otherwise, setup a
    # corresponding dictionary that returns a random float for each strain.
    if args.priority:
        priorities = read_priority_scores(args.priority)
    else:
        random_generator = np.random.default_rng(args.subsample_seed)
        priorities = defaultdict(random_generator.random)

    # Generate columns for grouping.
    grouping_metadata = enrich_metadata(
        metadata,
        group_by,
    )

    # Enrich with priorities.
    grouping_metadata['priority'] = [priorities[strain] for strain in grouping_metadata.index]

    pandas_groupby = grouping_metadata.groupby(list(group_by), group_keys=False)

    n_groups = len(pandas_groupby.groups)

    # Determine sequences per group.
    if args.sequences_per_group:
        sequences_per_group = args.sequences_per_group
    elif args.subsample_max_sequences:
        group_sizes = [len(strains) for strains in pandas_groupby.groups.values()]

        try:
            # Calculate sequences per group. If there are more groups than maximum
            # sequences requested, sequences per group will be a floating point
            # value and subsampling will be probabilistic.
            sequences_per_group, probabilistic_used = calculate_sequences_per_group(
                args.subsample_max_sequences,
                group_sizes,
                allow_probabilistic=args.probabilistic_sampling
            )
        except TooManyGroupsError as error:
            raise AugurError(str(error)) from error

        if (probabilistic_used):
            print(f"Sampling probabilistically at {sequences_per_group:0.4f} sequences per group, meaning it is possible to have more than the requested maximum of {args.subsample_max_sequences} sequences after filtering.")
        else:
            print(f"Sampling at {sequences_per_group} per group.")
    else:
        pass
        # FIXME: what to do when no subsampling is requested?

    group_size_limits = (size for size in get_group_size_limits(n_groups, sequences_per_group, random_seed=args.subsample_seed))

    def row_sampler(group):
        n = next(group_size_limits)
        return group.nlargest(n, 'priority')

    return {strain for strain in pandas_groupby.apply(row_sampler).index}

def enrich_metadata(metadata, group_by=None):
    """Enrich metadata with generated columns.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata to inspect for the given strains.
    group_by : list
        A list of metadata (or generated) columns to group records by.

    Returns
    -------
    metadata : pandas.DataFrame
        Metadata with generated columns.

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
    {'strain1': (2020, (2020, 1)), 'strain2': (2020, (2020, 2))}

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
    {'strain1': (2020, (2020, 1), 'unknown'), 'strain2': (2020, (2020, 2), 'unknown')}

    We can group metadata without any non-ID columns.

    >>> metadata = pd.DataFrame([{"strain": "strain1"}, {"strain": "strain2"}]).set_index("strain")
    >>> get_groups_for_subsampling(strains, metadata, group_by=('_dummy',))
    {'strain1': ('_dummy',), 'strain2': ('_dummy',)}
    """
    if len(metadata) == 0:
        return metadata

    if not group_by or group_by == ('_dummy',):
        metadata['_dummy'] = '_dummy'
        return metadata

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
            df_dates = metadata[METADATA_DATE_COLUMN].str.split('-', n=2, expand=True)
            df_dates = df_dates.set_axis(temp_date_cols[:len(df_dates.columns)], axis=1)
            missing_date_cols = set(temp_date_cols) - set(df_dates.columns)
            for col in missing_date_cols:
                df_dates[col] = pd.NA
            for col in temp_date_cols:
                df_dates[col] = pd.to_numeric(df_dates[col], errors='coerce').astype(pd.Int64Dtype())

            # Extend metadata with generated date columns
            # Drop the date column since it should not be used for grouping.
            metadata = pd.concat([metadata.drop(METADATA_DATE_COLUMN, axis=1), df_dates], axis=1)

            # FIXME: I think this is useless - drop it in another commit
            # Check again if metadata is empty after dropping ambiguous dates.
            # if metadata.empty:
            #     return group_by_strain

            # Generate columns.
            if constants.DATE_YEAR_COLUMN in generated_columns_requested:
                metadata[constants.DATE_YEAR_COLUMN] = metadata[f'{temp_prefix}year']
            if constants.DATE_MONTH_COLUMN in generated_columns_requested:
                metadata[constants.DATE_MONTH_COLUMN] = list(zip(
                    metadata[f'{temp_prefix}year'],
                    metadata[f'{temp_prefix}month']
                ))
            if constants.DATE_WEEK_COLUMN in generated_columns_requested:
                # Note that week = (year, week) from the date.isocalendar().
                # Do not combine the raw year with the ISO week number alone,
                # since raw year ≠ ISO year.
                metadata[constants.DATE_WEEK_COLUMN] = metadata.apply(lambda row: get_iso_year_week(
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

    return metadata


def get_group_size_limits(number_of_groups: int, max_size, max_attempts = 100, random_seed = None):
    """Return a list of group size limits.

    When the maximum size is fractional, probabilistically sample the maximum
    size from a Poisson distribution. Make at least the given number of maximum
    attempts to create groups for which the sum of their maximum sizes is
    greater than zero.

    Parameters
    ----------
    number_of_groups : int
        The number of groups.
    max_size : int | float
        Maximum size of a group.
    max_attempts : int
        Maximum number of attempts for creating group sizes.
    random_seed
        Seed value for np.random.default_rng for reproducible randomness.
    """
    sizes = None
    total_max_size = 0
    attempts = 0

    # If max_size is not fractional, use it as the limit for all groups.
    if int(max_size) == max_size:
        return np.full(number_of_groups, max_size)

    # For small fractional maximum sizes, it is possible to randomly select
    # maximum sizes that all equal zero. When this happens, filtering
    # fails unexpectedly. We make multiple attempts to create sizes with
    # maximum sizes greater than zero for at least one group.
    while total_max_size == 0 and attempts < max_attempts:
        sizes = np.random.default_rng(random_seed).poisson(max_size, size=number_of_groups)
        total_max_size = sum(sizes)
        attempts += 1

    assert sizes is not None

    return sizes


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
            sequences_per_group = _calculate_fractional_sequences_per_group(
                target_max_value,
                group_sizes,
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


def _calculate_fractional_sequences_per_group(
        target_max_value: int,
        group_sizes: Collection[int]
) -> float:
    """Returns the fractional sequences per group for the given list of group
    sequences such that the total doesn't exceed the requested number of
    samples.

    Parameters
    ----------
    target_max_value : int
        the total number of sequences allowed across all groups
    group_sizes : Collection[int]
        the number of sequences in each group

    Returns
    -------
    float
        fractional maximum number of sequences allowed per group to meet the
        required maximum total sequences allowed

    Examples
    --------
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
    hi = float(target_max_value)

    while (hi / lo) > 1.1:
        mid = (lo + hi) / 2
        if _calculate_total_sequences(mid, group_sizes) <= target_max_value:
            lo = mid
        else:
            hi = mid

    return (lo + hi) / 2
