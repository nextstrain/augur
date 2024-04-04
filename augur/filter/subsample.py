import numpy as np
from typing import Collection, Iterable, List, Sequence, Set
from augur.errors import AugurError
from augur.filter.debug import add_debugging
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.print import print_err
from augur.io.sqlite3 import Sqlite3Database, sanitize_identifier
from augur.io.tabular_file import TabularFile
from . import constants
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


    if args.weights:
        assert args.subsample_max_sequences is not None

        group_by_columns = list(TabularFile(args.weights, delimiter='\t').columns)
        group_by_columns.remove('weight')

        valid_group_by_columns = get_valid_group_by_columns(metadata_columns, group_by_columns)
        create_grouping_table(valid_group_by_columns)
        create_group_size_limits_table_weighted(args.weights, args.subsample_max_sequences)

    else:
        valid_group_by_columns = []
        if args.group_by:
            valid_group_by_columns = get_valid_group_by_columns(metadata_columns, args.group_by)
        create_grouping_table(valid_group_by_columns)

        if args.subsample_max_sequences:
            if args.group_by:
                group_sizes = get_group_sizes(valid_group_by_columns)
            else:
                valid_group_by_columns = [constants.GROUP_BY_DUMMY_COLUMN]
                group_sizes = [_get_filtered_strains_count()]

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
            sequences_per_group = args.sequences_per_group

        create_group_size_limits_table(valid_group_by_columns, sequences_per_group, args.subsample_seed)

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


def create_grouping_table(group_by_columns: Iterable[str]):
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
                TRUE AS {constants.GROUP_BY_DUMMY_COLUMN}
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


def get_group_sizes(group_by_columns: Iterable[str]):
    """
    Returns
    -------
    list of int
        A list with the number of sequences in each group.
    """
    group_by_columns_sql = ','.join(sanitize_identifier(column) for column in group_by_columns)
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT
                COUNT(*) AS count,
                {group_by_columns_sql}
            FROM {constants.GROUPING_TABLE}
            GROUP BY {group_by_columns_sql}
        """)
        return [int(row['count']) for row in result]


def create_group_size_limits_table(group_by_columns: Sequence[str], sequences_per_group: float, random_seed):
    """Create a table for group size limits."""
    number_of_groups = _get_number_of_groups(group_by_columns)

    group_size_limits = get_group_size_limits(number_of_groups, sequences_per_group, random_seed=random_seed)
    assert len(group_size_limits) == number_of_groups

    # Create a function to iterate over group size limits using consecutive calls.
    # This allows the corresponding column in the table to be created with one
    # value per row.
    generator_index = -1
    def group_size_limit_iterator():
        """Return the next group size limit."""
        nonlocal generator_index
        generator_index += 1

        # Since group_size_limits is an array of numpy integers, int() is used to ensure it is stored properly in the database.
        return int(group_size_limits[generator_index])

    group_by_columns_sql = ','.join(sanitize_identifier(column) for column in group_by_columns)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        # Register SQLite3 user-defined function.
        db.connection.create_function(group_size_limit_iterator.__name__, 0, group_size_limit_iterator)

        db.connection.execute(f"""CREATE TABLE {constants.GROUP_SIZE_LIMITS_TABLE} AS
            SELECT
                {group_by_columns_sql},
                {group_size_limit_iterator.__name__}() AS {constants.GROUP_SIZE_LIMIT_COLUMN}
            FROM {constants.GROUPING_TABLE}
            GROUP BY {group_by_columns_sql}
        """)

        # Remove user-defined function.
        db.connection.create_function(group_size_limit_iterator.__name__, 0, None)


def create_group_size_limits_table_weighted(weights_file: str, total_size: int):
    import pandas as pd
    weights = pd.read_csv(weights_file, delimiter='\t')
    total_weights = weights['weight'].sum()
    weights['__augur_size'] = total_size * weights['weight'] / total_weights
    group_by_columns = weights.columns.to_list()
    group_by_columns.remove('weight')
    group_by_columns.remove('__augur_size')
    group_by_columns_sql = ','.join(sanitize_identifier(column) for column in group_by_columns)

    # FIXME: load weights into a table and JOIN instead
    def get_size(*values):
        weights2 = weights.copy()
        for column, value in zip(group_by_columns, values):
            # FIXME: work with generated date columns
            # if column == 'month':
            #     value = lambda x: return ()
            weights2 = weights2[weights2[column] == value]
        if len(weights2) == 0:
            # No weights defined for this combination. Use a weight of 0.
            return 0
        return weights2['__augur_size']

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        # Register SQLite3 user-defined function.
        db.connection.create_function(get_size.__name__, len(group_by_columns), get_size)

        db.connection.execute(f"""CREATE TABLE {constants.GROUP_SIZE_LIMITS_TABLE} AS
            SELECT
                {group_by_columns_sql},
                {get_size.__name__}({group_by_columns_sql}) AS {constants.GROUP_SIZE_LIMIT_COLUMN}
            FROM {constants.GROUPING_TABLE}
            GROUP BY {group_by_columns_sql}
        """)
    

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

def _get_number_of_groups(group_by_columns: Iterable[str]):
    """"""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT COUNT(*) AS count
            FROM (
                SELECT DISTINCT {','.join(sanitize_identifier(column) for column in group_by_columns)}
                FROM {constants.GROUPING_TABLE}
            )
        """)
        return int(result.fetchone()["count"])