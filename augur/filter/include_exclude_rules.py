import ast
import json
import re
import pandas as pd
import sqlite3
import sqlparse
from typing import Any, Callable, Dict, List, Optional, Tuple
from augur.errors import AugurError
from augur.index import ID_COLUMN as SEQUENCE_INDEX_ID_COLUMN
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.print import print_err
from augur.io.strains import read_strains
from augur.io.sqlite3 import Sqlite3Database, sanitize_identifier
from augur.io.vcf import is_vcf as filename_is_vcf
from . import constants
from .debug import add_debugging

try:
    # pandas ≥1.5.0 only
    PandasUndefinedVariableError = pd.errors.UndefinedVariableError  # type: ignore
except AttributeError:
    PandasUndefinedVariableError = pd.core.computation.ops.UndefinedVariableError  # type: ignore

# A SQL expression to represent strains that the filter applies to.
SqlExpression = str

# Named parameters used in the SQL expression.
SqlParameters = Dict[str, Any]

# The return value of a filter function.
FilterFunctionReturn = Tuple[SqlExpression, SqlParameters]

# A function to use for filtering. Parameters vary.
FilterFunction = Callable[..., FilterFunctionReturn]

# A dictionary of parameters to pass to the filter function.
FilterFunctionKwargs = Dict[str, Any]

# A pair of values to represent the filter function and parameters to use.
FilterOption = Tuple[FilterFunction, FilterFunctionKwargs]


def filter_by_exclude_all() -> FilterFunctionReturn:
    """Exclude all strains regardless of the given metadata content.

    This is a placeholder function that can be called as part of a generalized
    loop through all possible functions.
    """
    expression = 'True'
    parameters: SqlParameters = {}
    return expression, parameters


def filter_by_exclude(exclude_file) -> FilterFunctionReturn:
    """Exclude the given strains.

    Parameters
    ----------
    exclude_file : str
        Filename with strain names to exclude
    """
    excluded_strains = read_strains(exclude_file)
    return _filter_by_exclude_strains(excluded_strains)


def _filter_by_exclude_strains(strains) -> FilterFunctionReturn:
    """Exclude the given strains.

    Parameters
    ----------
    exclude_file : str
        Filename with strain names to exclude
    """
    quoted_strains = (f"'{strain}'" for strain in strains)
    expression = f"""
        {constants.ID_COLUMN} IN ({','.join(quoted_strains)})
    """
    parameters: SqlParameters = {}
    return expression, parameters


def parse_filter_query(query):
    """Parse an augur filter-style query and return the corresponding column,
    operator, and value for the query.

    Parameters
    ----------
    query : str
        augur filter-style query following the pattern of `"property=value"` or `"property!=value"`

    Returns
    -------
    str :
        Name of column to query
    str :
        Either '=' or '!=' to denote the operator for a SQLite3 WHERE expression.
    str :
        Value of column to query

    Examples
    --------
    >>> parse_filter_query("property=value")
    ('property', '=', 'value')
    >>> parse_filter_query("property!=value")
    ('property', '!=', 'value')

    """
    column, value = re.split(r'!?=', query)
    op = '='
    if "!=" in query:
        op = '!='

    return column, op, value


def filter_by_exclude_where(exclude_where) -> FilterFunctionReturn:
    """Exclude all strains that match the given exclusion query.

    Unlike pandas query syntax, exclusion queries should follow the pattern of
    `"property=value"` or `"property!=value"`. Additionally, this filter treats
    all values like lowercase strings, so we convert all values to strings first
    and then lowercase them before testing the given query.

    Parameters
    ----------
    exclude_where : str
        Filter query used to exclude strains
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_columns = db.columns(constants.METADATA_TABLE)
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)

    column, op, value = parse_filter_query(exclude_where)

    if column in metadata_columns:
        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {sanitize_identifier(metadata_id_column)}
                FROM {constants.METADATA_TABLE}
                WHERE lower({constants.METADATA_TABLE}.{sanitize_identifier(column)}) {op} lower(:value)
            )
        """
        parameters = {'value': value}
    else:
        # Skip the filter, if the requested column does not exist.
        expression = 'False'
        parameters = {}

    return expression, parameters


def filter_by_sqlite_query(query: str, column_types: Optional[Dict[str, str]] = None) -> FilterFunctionReturn:
    """Filter by any valid SQLite expression on the metadata.

    Strains that do *not* match the query will be excluded.

    Parameters
    ----------
    query : str
        SQL expression used to exclude strains
    column_types : str
        Dict mapping of data type
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)
        metadata_columns = set(db.columns(constants.METADATA_TABLE))

    if column_types is None:
        column_types = {}

    # Set columns for type conversion.
    variables = extract_potential_sqlite_variables(query)
    if variables is not None:
        columns = variables.intersection(metadata_columns)
    else:
        # Column extraction failed. Apply type conversion to all columns.
        columns = metadata_columns

    # If a type is not explicitly provided, try converting the column to numeric.
    # This should cover most use cases, since one common problem is that the
    # built-in data type inference when loading the DataFrame does not
    # support nullable numeric columns, so numeric comparisons won't work on
    # those columns. pd.to_numeric does proper conversion on those columns,
    # and will not make any changes to columns with other values.
    for column in columns:
        column_types.setdefault(column, 'numeric')

    # FIXME: Apply column_types.
    # It's not easy to change the type on the table schema.¹
    # Maybe using CAST? But that always takes place even if the conversion is lossy
    # and irreversible (i.e. no error handling options like pd.to_numeric).
    # ¹ <https://www.sqlite.org/lang_altertable.html#making_other_kinds_of_table_schema_changes>
    # ² <https://www.sqlite.org/lang_expr.html#castexpr>

    expression = f"""
        {constants.ID_COLUMN} IN (
            SELECT {sanitize_identifier(metadata_id_column)}
            FROM {constants.METADATA_TABLE}
            WHERE NOT ({query})
        )
    """
    parameters: SqlParameters = {}
    return expression, parameters


def filter_by_pandas_query(query: str, chunksize: int, column_types: Optional[Dict[str, str]] = None) -> FilterFunctionReturn:
    """Filter by a Pandas expression on the metadata.

    Note that this is inefficient compared to native SQLite queries, and is in place
    for backwards compatibility.

    Parameters
    ----------
    query : str
        Pandas query string used on a DataFrame representation of the metadata.
    column_types : str
        Dict mapping of data type
    chunksize : int
        Maximum number of metadata records to read into memory at a time.
        Increasing this number can speed up filtering at the cost of more memory used.
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)
        metadata_columns = set(db.columns(constants.METADATA_TABLE))

        # TODO: select only the columns used in the query
        metadata_chunks = pd.read_sql_query(f"""
            SELECT *
            FROM {constants.METADATA_TABLE}
        """, db.connection, chunksize=chunksize)

    if column_types is None:
        column_types = {}

    # Set columns for type conversion.
    variables = extract_pandas_query_variables(query)
    if variables is not None:
        columns = variables.intersection(metadata_columns)
    else:
        # Column extraction failed. Apply type conversion to all columns.
        columns = metadata_columns

    # If a type is not explicitly provided, try automatic conversion.
    for column in columns:
        column_types.setdefault(column, 'auto')

    excluded_strains = []

    for metadata_chunk in metadata_chunks:
        # Convert data types before applying the query.
        # NOTE: This can behave differently between different chunks of metadata,
        # but it's the best we can do.
        for column, dtype in column_types.items():
            if dtype == 'auto':
                # Try numeric conversion followed by boolean conversion.
                try:
                    # pd.to_numeric supports nullable numeric columns unlike pd.read_csv's
                    # built-in data type inference.
                    metadata_chunk[column] = pd.to_numeric(metadata_chunk[column], errors='raise')
                except:
                    try:
                        metadata_chunk[column] = metadata_chunk[column].map(_string_to_boolean)
                    except ValueError:
                        # If both conversions fail, column values are preserved as strings.
                        pass

            elif dtype == 'int':
                try:
                    metadata_chunk[column] = pd.to_numeric(metadata_chunk[column], errors='raise', downcast='integer')
                except ValueError as e:
                    raise AugurError(f"Failed to convert value in column {column!r} to int. {e}")
            elif dtype == 'float':
                try:
                    metadata_chunk[column] = pd.to_numeric(metadata_chunk[column], errors='raise', downcast='float')
                except ValueError as e:
                    raise AugurError(f"Failed to convert value in column {column!r} to float. {e}")
            elif dtype == 'bool':
                try:
                    metadata_chunk[column] = metadata_chunk[column].map(_string_to_boolean)
                except ValueError as e:
                    raise AugurError(f"Failed to convert value in column {column!r} to bool. {e}")
            elif dtype == 'str':
                metadata_chunk[column] = metadata_chunk[column].astype('str', errors='ignore')

        try:
            matches = metadata_chunk.query(query).index
        except Exception as e:
            if isinstance(e, PandasUndefinedVariableError):
                raise AugurError(f"Query contains a column that does not exist in metadata.") from e
            raise AugurError(f"Internal Pandas error when applying query:\n\t{e}\nEnsure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.") from e

        # Exclude strains that do not match the query.
        excluded_strains.extend(
            metadata_chunk.drop(matches)[metadata_id_column].values
        )

    return _filter_by_exclude_strains(excluded_strains)


def _string_to_boolean(s: str):
    """Convert a string to an optional boolean value.

    Raises ValueError if it cannot be converted.
    """
    if s.lower() == 'true':
        return True
    elif s.lower() == 'false':
        return False
    elif s == '':
        return None

    raise ValueError(f"Unable to convert {s!r} to a boolean value.")


def filter_by_ambiguous_date(date_column, ambiguity) -> FilterFunctionReturn:
    """Filter where values in the given date column have a given level of ambiguity.

    Determine ambiguity hierarchically such that, for example, an ambiguous
    month implicates an ambiguous day even when day information is available.

    Parameters
    ----------
    date_column : str
        The date column is already parsed beforehand. However, this is used to
        verify that the column exists and it is still beneficial to report in
        the output log as a kwarg.
    ambiguity : str
        Level of date ambiguity to filter by
    """

    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_has_date_column = date_column in db.columns(constants.METADATA_TABLE)

    if not metadata_has_date_column:
        expression = 'False'
    elif ambiguity == 'year':
        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {constants.ID_COLUMN}
                FROM {constants.DATE_TABLE}
                WHERE {constants.DATE_YEAR_COLUMN} IS NULL
            )
        """
    elif ambiguity == 'month':
        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {constants.ID_COLUMN}
                FROM {constants.DATE_TABLE}
                WHERE (
                    {constants.DATE_MONTH_COLUMN} IS NULL OR
                    {constants.DATE_YEAR_COLUMN}  IS NULL
                )
            )
        """
    else:
        assert ambiguity == 'day' or ambiguity == 'any'
        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {constants.ID_COLUMN}
                FROM {constants.DATE_TABLE}
                WHERE (
                    {constants.DATE_DAY_COLUMN}   IS NULL OR
                    {constants.DATE_MONTH_COLUMN} IS NULL OR
                    {constants.DATE_YEAR_COLUMN}  IS NULL
                )
            )
        """
    parameters: SqlParameters = {}
    return expression, parameters


def skip_group_by_with_ambiguous_year(date_column) -> FilterFunctionReturn:
    """Alias to filter_by_ambiguous_date for year. This is to have a named function available for the filter reason."""
    return filter_by_ambiguous_date(date_column, ambiguity="year")


def skip_group_by_with_ambiguous_month(date_column) -> FilterFunctionReturn:
    """Alias to filter_by_ambiguous_date for month. This is to have a named function available for the filter reason."""
    return filter_by_ambiguous_date(date_column, ambiguity="month")


def skip_group_by_with_ambiguous_day(date_column) -> FilterFunctionReturn:
    """Alias to filter_by_ambiguous_date for day. This is to have a named function available for the filter reason."""
    return filter_by_ambiguous_date(date_column, ambiguity="day")


def filter_by_min_date(date_column, min_date) -> FilterFunctionReturn:
    """Filter by minimum date.

    Parameters
    ----------
    date_column
        The date column is already parsed beforehand. However, this is used to
        verify that the column exists and it is still beneficial to report in
        the output log as a kwarg.
    min_date : float
        Minimum date
    """

    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_has_date_column = date_column in db.columns(constants.METADATA_TABLE)

    # Skip this filter if the date column does not exist.
    if not metadata_has_date_column:
        expression = 'False'
        parameters: SqlParameters = {}
    else:
        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {constants.ID_COLUMN}
                FROM {constants.DATE_TABLE}
                WHERE {constants.NUMERIC_DATE_MAX_COLUMN} < :min_date OR {constants.NUMERIC_DATE_MIN_COLUMN} IS NULL
            )
        """
        parameters = {'min_date': min_date}

    return expression, parameters


def filter_by_max_date(date_column, max_date) -> FilterFunctionReturn:
    """Filter by maximum date.

    Parameters
    ----------
    date_column : str
        The date column is already parsed beforehand. However, this is used to
        verify that the column exists and it is still beneficial to report in
        the output log as a kwarg.
    max_date : float
        Maximum date
    """

    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_has_date_column = date_column in db.columns(constants.METADATA_TABLE)

    # Skip this filter if the date column does not exist.
    if not metadata_has_date_column:
        expression = 'False'
        parameters: SqlParameters = {}
    else:
        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {constants.ID_COLUMN}
                FROM {constants.DATE_TABLE}
                WHERE {constants.NUMERIC_DATE_MIN_COLUMN} > :max_date OR {constants.NUMERIC_DATE_MAX_COLUMN} IS NULL
            )
        """
        parameters = {'max_date': max_date}

    return expression, parameters


def filter_by_sequence_index() -> FilterFunctionReturn:
    """Filter by presence of corresponding entries in a given sequence
    index. This filter effectively intersects the strain ids in the metadata and
    sequence index.
    """
    expression = f"""
        {constants.ID_COLUMN} NOT IN (
            SELECT {SEQUENCE_INDEX_ID_COLUMN}
            FROM {constants.SEQUENCE_INDEX_TABLE}
        )
    """
    parameters: SqlParameters = {}
    return expression, parameters


def filter_by_sequence_length(min_length) -> FilterFunctionReturn:
    """Filter by sequence length from a given sequence index.

    Parameters
    ----------
    min_length : int
        Minimum number of standard nucleotide characters (A, C, G, or T) in each sequence
    """
    expression = f"""
        {constants.ID_COLUMN} IN (
            SELECT {SEQUENCE_INDEX_ID_COLUMN}
            FROM {constants.SEQUENCE_INDEX_TABLE}
            WHERE A + C + G + T < :min_length
        )
    """
    parameters = {'min_length': min_length}
    return expression, parameters


def filter_by_non_nucleotide() -> FilterFunctionReturn:
    """Filter for strains with invalid nucleotide content.
    """
    expression = f"""
        {constants.ID_COLUMN} IN (
            SELECT {SEQUENCE_INDEX_ID_COLUMN}
            FROM {constants.SEQUENCE_INDEX_TABLE}
            WHERE invalid_nucleotides != 0
        )
    """
    parameters: SqlParameters = {}
    return expression, parameters


def force_include_strains(include_file) -> FilterFunctionReturn:
    """Include strains in the given text file.

    Parameters
    ----------
    include_file : str
        Filename with strain names to include
    """
    strains = read_strains(include_file)
    quoted_strains = (f"'{strain}'" for strain in strains)
    expression = f"""
        {constants.ID_COLUMN} IN ({','.join(quoted_strains)})
    """
    parameters: SqlParameters = {}
    return expression, parameters


def force_include_where(include_where) -> FilterFunctionReturn:
    """Include all strains that match the given query.

    Unlike pandas query syntax, inclusion queries should follow the pattern of
    `"property=value"` or `"property!=value"`. Additionally, this filter treats
    all values like lowercase strings, so we convert all values to strings first
    and then lowercase them before testing the given query.

    Parameters
    ----------
    include_where : str
        Filter query used to include strains
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_columns = db.columns(constants.METADATA_TABLE)
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)

    column, op, value = parse_filter_query(include_where)

    if column in metadata_columns:

        expression = f"""
            {constants.ID_COLUMN} IN (
                SELECT {sanitize_identifier(metadata_id_column)}
                FROM {constants.METADATA_TABLE}
                WHERE {constants.METADATA_TABLE}.{sanitize_identifier(column)} {op} :value
            )
        """
        parameters = {'value': value}
    else:
        # Skip the inclusion filter if the requested column does not exist.
        expression = 'False'
        parameters = {}

    return expression, parameters


def construct_filters(args) -> Tuple[List[FilterOption], List[FilterOption]]:
    """Construct lists of filters and inclusion criteria based on user-provided
    arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments provided by the user.
    """
    exclude_by: List[FilterOption] = []
    include_by: List[FilterOption] = []

    # Force include sequences specified in file(s).
    if args.include:
        # Collect the union of all given strains to include.
        for include_file in args.include:
            include_by.append((
                force_include_strains,
                {
                    "include_file": include_file,
                }
            ))

    # Add sequences with particular metadata attributes.
    if args.include_where:
        for include_where in args.include_where:
            include_by.append((
                force_include_where,
                {
                    "include_where": include_where,
                }
            ))

    # Exclude all strains by default.
    if args.exclude_all:
        exclude_by.append((filter_by_exclude_all, {}))

    # Filter by sequence index.
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        if constants.SEQUENCE_INDEX_TABLE in db.tables():
            exclude_by.append((filter_by_sequence_index, {}))

    # Remove strains explicitly excluded by name.
    if args.exclude:
        for exclude_file in args.exclude:
            exclude_by.append((
                filter_by_exclude,
                {
                    "exclude_file": exclude_file,
                }
            ))

    # Exclude strain my metadata field like 'host=camel'.
    if args.exclude_where:
        for exclude_where in args.exclude_where:
            exclude_by.append((
                filter_by_exclude_where,
                {"exclude_where": exclude_where}
            ))

    # Exclude strains by metadata.
    if args.query_pandas:
        kwargs = {
            "query": args.query_pandas,
            "chunksize": args.metadata_chunk_size,
        }
        if args.query_columns:
            kwargs["column_types"] = {column: dtype for column, dtype in args.query_columns}

        exclude_by.append((
            filter_by_pandas_query,
            kwargs
        ))
    if args.query_sqlite:
        kwargs = {"query": args.query_sqlite}
        if args.query_columns:
            kwargs["column_types"] = {column: dtype for column, dtype in args.query_columns}

        exclude_by.append((
            filter_by_sqlite_query,
            kwargs
        ))

    # Filter by ambiguous dates.
    if args.exclude_ambiguous_dates_by:
        exclude_by.append((
            filter_by_ambiguous_date,
            {
                "date_column": METADATA_DATE_COLUMN,
                "ambiguity": args.exclude_ambiguous_dates_by,
            }
        ))

    # Filter by min/max date.
    if args.min_date:
        exclude_by.append((
            filter_by_min_date,
            {
                "min_date": args.min_date,
                "date_column": METADATA_DATE_COLUMN,
            }
        ))
    if args.max_date:
        exclude_by.append((
            filter_by_max_date,
            {
                "max_date": args.max_date,
                "date_column": METADATA_DATE_COLUMN,
            }
        ))

    # Filter by sequence length.
    if args.min_length:
        # Skip VCF files and warn the user that the min length filter does not
        # make sense for VCFs.
        is_vcf = filename_is_vcf(args.sequences)

        if is_vcf: #doesn't make sense for VCF, ignore.
            print_err("WARNING: Cannot use min_length for VCF files. Ignoring...")
        else:
            exclude_by.append((
                filter_by_sequence_length,
                {
                    "min_length": args.min_length,
                }
            ))

    # Exclude sequences with non-nucleotide characters.
    if args.non_nucleotide:
        exclude_by.append((filter_by_non_nucleotide, {}))

    if args.group_by:
        # The order in which these are applied later should be broad → specific
        # ambiguity (e.g. year then month), otherwise broad ambiguity will be
        # captured by specific ambiguity.
        if {constants.DATE_YEAR_COLUMN, constants.DATE_MONTH_COLUMN, constants.DATE_WEEK_COLUMN} & set(args.group_by):
            exclude_by.append((
                skip_group_by_with_ambiguous_year,
                {"date_column": METADATA_DATE_COLUMN}
            ))
        if {constants.DATE_MONTH_COLUMN, constants.DATE_WEEK_COLUMN} & set(args.group_by):
            exclude_by.append((
                skip_group_by_with_ambiguous_month,
                {"date_column": METADATA_DATE_COLUMN}
            ))
        if constants.DATE_WEEK_COLUMN in args.group_by:
            exclude_by.append((
                skip_group_by_with_ambiguous_day,
                {"date_column": METADATA_DATE_COLUMN}
            ))

    return exclude_by, include_by


@add_debugging
def apply_filters(exclude_by: List[FilterOption], include_by: List[FilterOption]):
    """Apply exclusion and force-inclusion rules to filter strains from the metadata."""
    init_filter_reason_table()
    apply_exclusions(exclude_by)
    apply_force_inclusions(include_by)


def init_filter_reason_table():
    """Initialize the filter reason table with all strains as not being excluded nor force-included."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.connection.execute(f"""CREATE TABLE {constants.FILTER_REASON_TABLE} AS
            SELECT
                {sanitize_identifier(metadata_id_column)} AS {constants.ID_COLUMN},
                FALSE AS {constants.EXCLUDE_COLUMN},
                FALSE AS {constants.INCLUDE_COLUMN},
                NULL AS {constants.FILTER_REASON_COLUMN},
                NULL AS {constants.FILTER_REASON_KWARGS_COLUMN}
            FROM {constants.METADATA_TABLE}
        """)
        db.create_primary_index(constants.FILTER_REASON_TABLE, constants.ID_COLUMN)


def apply_exclusions(exclude_by: List[FilterOption]):
    """Update the filter reason table given the outcome of exclusion filters."""

    # Reversed so that earlier entries in the original list will be have higher precedence.
    # This is because later evaluations will overwrite any previously applied filter reasons.
    for exclude_function, kwargs in reversed(exclude_by):
        where_expression, where_parameters = None, None

        # Note: Consider using JOIN instead of subqueries if performance issues arise¹.
        # ¹ https://stackoverflow.com/q/3856164
        where_expression, where_parameters = exclude_function(**kwargs)

        assert where_expression is not None
        assert where_parameters is not None

        sql = f"""
            UPDATE {constants.FILTER_REASON_TABLE}
            SET
                {constants.EXCLUDE_COLUMN} = TRUE,
                {constants.FILTER_REASON_COLUMN} = :filter_reason,
                {constants.FILTER_REASON_KWARGS_COLUMN} = :filter_reason_kwargs
            WHERE {where_expression}
        """

        sql_parameters = {
            'filter_reason': exclude_function.__name__,
            'filter_reason_kwargs': filter_kwargs_to_str(kwargs)
        }

        # Add parameters returned from the filter function.
        sql_parameters = {**sql_parameters, **where_parameters}

        with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
            try:
                db.connection.execute(sql, sql_parameters)
            except Exception as e:
                if exclude_function is filter_by_sqlite_query:
                    if isinstance(e, sqlite3.OperationalError):
                        if "no such column" in str(e):
                            raise AugurError(f"Query contains a column that does not exist in metadata.") from e
                        raise AugurError(f"Error when applying query. Ensure the syntax is valid per <https://www.sqlite.org/lang_expr.html>.") from e


def apply_force_inclusions(include_by: List[FilterOption]):
    """Update the filter reason table with force-inclusion rules."""
    for include_function, kwargs in include_by:
        where_expression, where_parameters = include_function(**kwargs)
        sql = f"""
            UPDATE {constants.FILTER_REASON_TABLE}
            SET
                {constants.INCLUDE_COLUMN} = TRUE,
                {constants.FILTER_REASON_COLUMN} = :filter_reason,
                {constants.FILTER_REASON_KWARGS_COLUMN} = :filter_reason_kwargs
            WHERE {where_expression}
        """

        sql_parameters = {
            'filter_reason': include_function.__name__,
            'filter_reason_kwargs': filter_kwargs_to_str(kwargs)
        }

        # Add parameters returned from the filter function.
        sql_parameters = {**sql_parameters, **where_parameters}

        with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
            db.connection.execute(sql, sql_parameters)


def filter_kwargs_to_str(kwargs: FilterFunctionKwargs):
    """Convert a dictionary of kwargs to a JSON string for downstream reporting.

    This structured string can be converted back into a Python data structure
    later for more sophisticated reporting by specific kwargs.

    This function excludes data types from arguments like pandas DataFrames and
    also converts floating point numbers to a fixed precision for better
    readability and reproducibility.

    Parameters
    ----------
    kwargs : dict
        Dictionary of kwargs passed to a given filter function.

    Returns
    -------
    str :
        String representation of the kwargs for reporting.

    Examples
    --------
    >>> from augur.dates import numeric_date
    >>> from augur.filter.include_exclude_rules import filter_by_sequence_length, filter_by_min_date
    >>> exclude_by = [(filter_by_sequence_length, {"min_length": 27000})]
    >>> filter_kwargs_to_str(exclude_by[0][1])
    '[["min_length", 27000]]'
    >>> exclude_by = [(filter_by_min_date, {"date_column": "date", "min_date": numeric_date("2020-03-01")})]
    >>> filter_kwargs_to_str(exclude_by[0][1])
    '[["date_column", "date"], ["min_date", 2020.17]]'

    """
    # Sort keys prior to processing to guarantee the same output order
    # regardless of the input order.
    sorted_keys = sorted(kwargs.keys())

    kwarg_list = []
    for key in sorted_keys:
        value = kwargs[key]

        # Handle special cases for data types that we want to represent
        # differently from their defaults.
        if isinstance(value, float):
            value = round(value, 2)

        # Don't include chunksize since it does not affect end results.
        if key == 'chunksize':
            continue

        kwarg_list.append((key, value))

    return json.dumps(kwarg_list)


def extract_pandas_query_variables(pandas_query: str):
    """Try extracting all variable names used in a pandas query string.

    If successful, return the variable names as a set. Otherwise, nothing is returned.

    Examples
    --------
    >>> extract_pandas_query_variables("var1 == 'value'")
    {'var1'}
    >>> sorted(extract_pandas_query_variables("var1 == 'value' & var2 == 10"))
    ['var1', 'var2']
    >>> extract_pandas_query_variables("var1.str.startswith('prefix')")
    {'var1'}
    >>> extract_pandas_query_variables("this query is invalid")
    """
    # Since Pandas' query grammar should be a subset of Python's, which uses the
    # ast stdlib under the hood, we can try to parse queries with that as well.
    # Errors may arise from invalid query syntax or any Pandas syntax not
    # covered by Python (unlikely, but I'm not sure). In those cases, don't
    # return anything.
    try:
        return set(node.id
                   for node in ast.walk(ast.parse(pandas_query))
                   if isinstance(node, ast.Name))
    except:
        return None


def extract_potential_sqlite_variables(sqlite_expression: str):
    """Try extracting all variable names used in a SQLite expression.

    If successful, return the variable names as a set. Otherwise, nothing is returned.

    Examples
    --------
    >>> extract_potential_sqlite_variables("var1 = 'value'")
    {'var1'}
    >>> sorted(extract_potential_sqlite_variables("var1 = 'value' AND var2 = 10"))
    ['var1', 'var2']
    >>> extract_potential_sqlite_variables("var1 LIKE 'prefix%'")
    {'var1'}
    >>> sorted(extract_potential_sqlite_variables("this query is invalid"))
    ['invalid', 'this query']
    """
    # This seems to be more difficult than Pandas query parsing.
    # <https://stackoverflow.com/q/35624662>
    try:
        query = f"SELECT * FROM table WHERE {sqlite_expression}"
        where = [x for x in sqlparse.parse(query)[0] if isinstance(x, sqlparse.sql.Where)][0]
        variables = set(_get_identifiers(where)) or None
        return variables
    except:
        return None


def _get_identifiers(token: sqlparse.sql.Token):
    """Yield identifiers from a token's children.
    
    Inspired by ast.walk.
    """
    from collections import deque
    todo = deque([token])
    while todo:
        node = todo.popleft()

        # Limit to comparisons to avoid false positives.
        # I chose not to use this because it also comes with false negatives.
        #
        # if isinstance(node, sqlparse.sql.Comparison):
        #     if isinstance(node.left, sqlparse.sql.Identifier):
        #         yield str(node.left)
        #     elif hasattr(node.left, 'tokens'):
        #         todo.extend(node.left.tokens)
        #     if isinstance(node.right, sqlparse.sql.Identifier):
        #         yield str(node.right)
        #     elif hasattr(node.right, 'tokens'):
        #         todo.extend(node.right.tokens)

        if isinstance(node, sqlparse.sql.Identifier):
            yield str(node)
        elif hasattr(node, 'tokens'):
            todo.extend(node.tokens)
