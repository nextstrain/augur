from collections import defaultdict
import itertools
import numpy as np
import pandas as pd
from textwrap import dedent
from typing import Collection, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from augur.errors import AugurError
from augur.filter.debug import add_debugging
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.print import print_err
from augur.io.sqlite3 import Sqlite3Database, sanitize_identifier
from . import constants
from .weights_file import WEIGHTS_COLUMN, COLUMN_VALUE_FOR_DEFAULT_WEIGHT, get_default_weight, get_weighted_columns, read_weights_file

Group = Tuple[str, ...]
"""Combination of grouping column values in tuple form."""
from .io import import_priorities_table


def get_valid_group_by_columns(metadata_columns: Set[str], group_by: List[str]):
    """Perform validation on requested group-by columns and return the valid subset.

    Parameters
    ----------
    metadata_columns
        All column names in metadata.
    group_by
        A list of metadata (or generated) columns to group records by.

    Returns
    -------
    list of str:
        Valid group-by columns.
    """
    # Create a set copy for faster existence checks.
    group_by_set = set(group_by)

    generated_columns_requested = constants.GROUP_BY_GENERATED_COLUMNS & group_by_set

    # If we could not find any requested categories, we cannot complete subsampling.
    if METADATA_DATE_COLUMN not in metadata_columns and group_by_set <= constants.GROUP_BY_GENERATED_COLUMNS:
        raise AugurError(f"The specified group-by categories ({group_by}) were not found. Note that using any of {sorted(constants.GROUP_BY_GENERATED_COLUMNS)} requires a column called {METADATA_DATE_COLUMN!r}.")
    if not group_by_set & (set(metadata_columns) | constants.GROUP_BY_GENERATED_COLUMNS):
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

        if METADATA_DATE_COLUMN not in metadata_columns:
            print_err(f"WARNING: A {METADATA_DATE_COLUMN!r} column could not be found to group-by {sorted(generated_columns_requested)}.")
            print_err(f"Filtering by group may behave differently than expected!")

    unknown_groups = group_by_set - metadata_columns - constants.GROUP_BY_GENERATED_COLUMNS
    if unknown_groups:
        print_err(f"WARNING: Some of the specified group-by categories couldn't be found: {', '.join(unknown_groups)}")
        print_err("Filtering by group may behave differently than expected!")
        for group in unknown_groups:
            group_by.remove(group)
    return group_by


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


# FIXME: read weighs file into sql table?

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
            sequences = 'sequence' if row[TARGET_SIZE_COLUMN] == 1 else 'sequences'
            are = 'is' if row[INPUT_SIZE_COLUMN] == 1 else 'are'
            group = list(f'{col}={value!r}' for col, value in row[group_by].items())
            print_err(f"WARNING: Targeted {row[TARGET_SIZE_COLUMN]} {sequences} for group {group} but only {row[INPUT_SIZE_COLUMN]} {are} available.")

    if output_sizes_file:
        # FIXME: make the order of rows deterministic
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
    lists = [list(values_for_unweighted_columns[column]) for column in unweighted_columns]
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
        unit = "group" if count == 1 else "groups"
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
    columns = 'column' if len(unweighted_columns) == 1 else 'columns'
    those = 'that' if len(unweighted_columns) == 1 else 'those'
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
            The default weight of {default_weight!r} will be used for all groups defined by those values."""))

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


@add_debugging
def apply_subsampling(args):
    """Apply subsampling to update the filter reason table.

    We handle the following major use cases:

    1. group by and sequences per group defined -> use the given values by the
    user to identify the highest priority records from each group.

    2. group by and maximum sequences defined -> count the group sizes, calculate the
    sequences per group that satisfies the requested maximum, and select that many sequences per group.

    3. group by not defined but maximum sequences defined -> use a "dummy"
    group such that we select at most the requested maximum number of
    sequences.
    """

    # Each strain has a score to determine priority during subsampling.
    # When no priorities are provided, they will be randomly generated.
    create_priorities_table(args)

    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_columns = set(db.columns(constants.METADATA_TABLE))

    # FIXME: optimize conditions

    valid_group_by_columns = []
    if args.group_by:
        valid_group_by_columns = get_valid_group_by_columns(metadata_columns, args.group_by)
        
    create_grouping_table(valid_group_by_columns, metadata_columns)

    if not args.group_by:
        valid_group_by_columns = [constants.GROUP_BY_DUMMY_COLUMN]
        target_group_sizes = {(constants.GROUP_BY_DUMMY_VALUE, ): _get_filtered_strains_count()}

    records_per_group = get_records_per_group(valid_group_by_columns)

    if args.subsample_max_sequences:
            if args.group_by_weights:
                print_err(f"Sampling with weights defined by {args.group_by_weights}.")
                target_group_sizes = get_weighted_group_sizes(
                    records_per_group,
                    args.group_by,
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
                    target_group_sizes = get_probabilistic_group_sizes(
                        records_per_group.keys(),
                        sequences_per_group,
                        random_seed=args.subsample_seed,
                    )
                else:
                    print_err(f"Sampling at {sequences_per_group} per group.")
                    assert type(sequences_per_group) is int
                    target_group_sizes = {group: sequences_per_group for group in records_per_group.keys()}
    else:
        assert args.sequences_per_group
        target_group_sizes = {group: args.sequences_per_group for group in records_per_group.keys()}

    create_group_size_limits_table(valid_group_by_columns, target_group_sizes)
    update_filter_reason_table(valid_group_by_columns)


def create_priorities_table(args):
    """Import or generate the priorities table."""
    if args.priority:
        import_priorities_table(args.priority)
    else:
        generate_priorities_table(args.subsample_seed)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.create_primary_index(constants.PRIORITIES_TABLE, constants.ID_COLUMN)


def generate_priorities_table(random_seed: int = None):
    """Generate a priorities table with random priorities.

    It is not possible to seed the SQLite built-in RANDOM(). As an alternative,
    use a Python function registered as a user-defined function.

    The generated priorities are random floats in the half-open interval [0.0, 1.0).
    """
    rng = np.random.default_rng(random_seed)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        # Register SQLite3 user-defined function.
        db.connection.create_function(rng.random.__name__, 0, rng.random)

        db.connection.execute(f"""CREATE TABLE {constants.PRIORITIES_TABLE} AS
            SELECT
                {constants.ID_COLUMN},
                {rng.random.__name__}() AS {constants.PRIORITY_COLUMN}
            FROM {constants.FILTER_REASON_TABLE}
            WHERE NOT {constants.EXCLUDE_COLUMN} OR {constants.INCLUDE_COLUMN}
        """)

        # Remove user-defined function.
        db.connection.create_function(rng.random.__name__, 0, None)


def create_grouping_table(group_by_columns: Iterable[str], metadata_columns: Set[str]):
    """Create a table with columns for grouping."""

    # For both of these, start with an empty string in case it isn't needed.
    generated_group_by_columns_sql = ''
    metadata_group_by_columns_sql = ''

    if group_by_columns:
        group_by_columns_set = set(group_by_columns)

        generated_group_by_columns = constants.GROUP_BY_GENERATED_COLUMNS & group_by_columns_set

        if generated_group_by_columns:
            generated_group_by_columns_sql = (
                # Prefix columns with the table alias defined in the SQL query further down.
                ','.join(f'd.{column}' for column in generated_group_by_columns)
                # Add an extra comma for valid SQL.
                + ',')

        metadata_group_by_columns = group_by_columns_set - generated_group_by_columns

        if metadata_group_by_columns:
            metadata_group_by_columns_sql = (
                # Prefix columns with the table alias defined in the SQL query further down.
                ','.join(f'm.{sanitize_identifier(column)}' for column in metadata_group_by_columns)
                # Add an extra comma for valid SQL.
                + ',')

    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)

    # Create a new table with rows as filtered metadata, with the following columns:
    # - Metadata ID column
    # - Group-by columns
    # - Generated date columns
    # - Priority score column
    # - Placeholder column
    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.connection.execute(f"""CREATE TABLE {constants.GROUPING_TABLE} AS
            SELECT
                f.{constants.ID_COLUMN},
                {metadata_group_by_columns_sql}
                {generated_group_by_columns_sql}
                p.{constants.PRIORITY_COLUMN},
                {constants.GROUP_BY_DUMMY_VALUE} AS {constants.GROUP_BY_DUMMY_COLUMN}
            FROM {constants.FILTER_REASON_TABLE} AS f
            JOIN {constants.METADATA_TABLE} AS m
                ON (f.{constants.ID_COLUMN} = m.{sanitize_identifier(metadata_id_column)})
            JOIN {constants.DATE_TABLE} AS d
                USING ({constants.ID_COLUMN})
            LEFT OUTER JOIN {constants.PRIORITIES_TABLE} AS p
                USING ({constants.ID_COLUMN})
            WHERE
                NOT f.{constants.EXCLUDE_COLUMN} OR f.{constants.INCLUDE_COLUMN}
        """)
    # Note: The last JOIN is a LEFT OUTER JOIN since a default INNER JOIN would
    # drop strains without a priority.


def get_records_per_group(group_by_columns: Sequence[str]) -> Dict[Group, int]:
    group_by_columns_sql = ','.join(sanitize_identifier(column) for column in group_by_columns)
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT
                COUNT(*) AS count,
                {group_by_columns_sql}
            FROM {constants.GROUPING_TABLE}
            GROUP BY {group_by_columns_sql}
        """)
        return {tuple(row[c] for c in group_by_columns) : int(row['count']) for row in result}


# FIXME: use group_sizes instead of iterator
def create_group_size_limits_table(group_by_columns: Sequence[str], target_group_sizes: Dict[Group, float]):
    """Create a table for group size limits."""

    # Create a function to return target group size for a given group
    def group_size(*group):
        return target_group_sizes[group]

    group_by_columns_sql = ','.join(sanitize_identifier(column) for column in group_by_columns)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        # Register SQLite3 user-defined function.
        db.connection.create_function(group_size.__name__, -1, group_size)

        db.connection.execute(f"""CREATE TABLE {constants.GROUP_SIZE_LIMITS_TABLE} AS
            SELECT
                {group_by_columns_sql},
                {group_size.__name__}({group_by_columns_sql}) AS {constants.GROUP_SIZE_LIMIT_COLUMN}
            FROM {constants.GROUPING_TABLE}
            GROUP BY {group_by_columns_sql}
        """)

        # Remove user-defined function.
        db.connection.create_function(group_size.__name__, 0, None)


def update_filter_reason_table(group_by_columns: Iterable[str]):
    """Subsample filtered metadata and update the filter reason table."""
    group_by_columns_sql = ','.join(sanitize_identifier(column) for column in group_by_columns)

    # First, select the strain column, group-by columns, and a `group_rank`
    # variable from the grouping table. `group_rank` represents an incremental
    # number ordered by priority within each group (i.e. the highest priority
    # strain per group gets group_rank=0).
    strains_with_group_rank = f"""
        SELECT
            {constants.ID_COLUMN},
            {group_by_columns_sql},
            ROW_NUMBER() OVER (
                PARTITION BY {group_by_columns_sql}
                ORDER BY
                    (CASE WHEN {constants.PRIORITY_COLUMN} IS NULL THEN 1 ELSE 0 END),
                    CAST({constants.PRIORITY_COLUMN} AS REAL)
                    DESC
            ) AS group_rank
        FROM {constants.GROUPING_TABLE}
    """
    # Notes:
    # 1. Although the name is similar to --group-by, the GROUP BY clause does not
    #    apply here. That command is used for aggregation commands such as getting
    #    the sizes of each group, which is done elsewhere.
    # 2. To treat rows without priorities as lowest priority, `ORDER BY … NULLS LAST`
    #    would be ideal. However, that syntax is unsupported on SQLite <3.30.0¹ so
    #    `CASE … IS NULL …` is a more widely compatible equivalent.
    #    ¹ https://www.sqlite.org/changes.html

    # Combine the above with the group size limits table to select the highest
    # priority strains.
    query_for_subsampled_strains = f"""
        SELECT {constants.ID_COLUMN}
        FROM ({strains_with_group_rank})
        JOIN {constants.GROUP_SIZE_LIMITS_TABLE} USING ({group_by_columns_sql})
        WHERE group_rank <= {constants.GROUP_SIZE_LIMIT_COLUMN}
    """

    # Exclude strains that didn't pass subsampling.
    # Note that the exclude column was already considered when creating the
    # grouping table earlier. The condition here is only in place to not
    # overwrite any existing reason for exclusion.
    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.connection.execute(f"""
            UPDATE {constants.FILTER_REASON_TABLE}
            SET
                {constants.EXCLUDE_COLUMN} = TRUE,
                {constants.FILTER_REASON_COLUMN} = '{constants.SUBSAMPLE_FILTER_REASON}'
            WHERE (
                NOT {constants.EXCLUDE_COLUMN}
                AND {constants.ID_COLUMN} NOT IN ({query_for_subsampled_strains})
            )
        """)


def _get_filtered_strains_count():
    """Returns the number of metadata strains that pass all filter rules."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT COUNT(*) AS count
            FROM {constants.FILTER_REASON_TABLE}
            WHERE NOT {constants.EXCLUDE_COLUMN} OR {constants.INCLUDE_COLUMN}
        """)
        return int(result.fetchone()["count"])
