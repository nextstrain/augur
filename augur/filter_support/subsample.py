import numpy as np
import pandas as pd
from typing import Collection, List, Set

from augur.filter_support.exceptions import FilterException
from augur.io import print_err


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
            print_err(f"WARNING: {error}")
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
    subsample.TooManyGroupsError: Asked to provide at most 1 sequences, but there are 2 groups.
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


def get_sizes_per_group(df_groups:pd.DataFrame, size_col:str, max_size, max_attempts=100, random_seed=None):
    """Create a DataFrame of sizes per group for the given maximum size.

    When the maximum size is fractional, probabilistically sample the maximum
    size from a Poisson distribution. Make at least the given number of maximum
    attempts to create groups for which the sum of their maximum sizes is
    greater than zero.

    Parameters
    ----------
    df_groups : pd.DataFrame
        DataFrame with one row per group.
    size_col : str
        Name of new column for group size.
    max_size : int | float
        Maximum size of a group.
    max_attempts : int
        Maximum number of attempts for creating group sizes.
    random_seed
        Seed value for np.random.default_rng for reproducible randomness.

    Returns
    -------
    pd.DataFrame
        df_groups with an additional column (size_col) containing group size.
    """
    total_max_size = 0
    attempts = 0

    # For small fractional maximum sizes, it is possible to randomly select
    # maximum sizes that all equal zero. When this happens, filtering
    # fails unexpectedly. We make multiple attempts to create sizes with
    # maximum sizes greater than zero for at least one group.
    while total_max_size == 0 and attempts < max_attempts:
        if max_size < 1.0:
            df_groups[size_col] = np.random.default_rng(random_seed).poisson(max_size, size=len(df_groups))
        else:
            df_groups[size_col] = max_size
        total_max_size = sum(df_groups[size_col])
        attempts += 1

    # TODO: raise error if total_max_size is still 0
    return df_groups


def get_valid_group_by_cols(group_by_cols:List[str], metadata_cols: Set[str]):
    """Perform validation on requested group-by columns and return the valid subset.

    Parameters
    ----------
    group_by_cols
        Column names requested for grouping.
    metadata_cols
        All column names in metadata.

    Returns
    -------
    list(str):
        Valid group-by columns.
    """
    # TODO: change behavior to address https://github.com/nextstrain/augur/issues/754
    extracted_date_cols = {'year', 'month'}
    group_by_set = set(group_by_cols)
    if group_by_set <= extracted_date_cols and 'date' not in metadata_cols:
        # all requested group-by columns are extracted date columns, but the date column is missing
        raise FilterException(f"The specified group-by categories ({group_by_cols}) were not found. Note that using 'year' or 'year month' requires a column called 'date'.")
    if not group_by_set & (metadata_cols | extracted_date_cols):
        # none of the requested group-by columns are valid
        raise FilterException(f"The specified group-by categories ({group_by_cols}) were not found.")
    unknown_cols = list(group_by_set - metadata_cols - extracted_date_cols)
    if 'date' not in metadata_cols:
        if "year" in group_by_set:
            print_err("WARNING: A 'date' column could not be found to group-by year.")
            unknown_cols.append("year")
        if "month" in group_by_set:
            print_err("WARNING: A 'date' column could not be found to group-by month.")
            unknown_cols.append("month")
    if unknown_cols:
        # warn and skip unknown columns
        print_err(f"WARNING: Some of the specified group-by categories couldn't be found: {', '.join(unknown_cols)}")
        print_err("Filtering by group may behave differently than expected!")
        valid_group_by_cols = list(group_by_cols)  # copy to preserve input object
        for col in unknown_cols:
            valid_group_by_cols.remove(col)
        return valid_group_by_cols
    return group_by_cols
