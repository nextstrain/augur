import argparse
import re
from typing import List, Set, Tuple
import numpy as np
import pandas as pd
import sqlite3
from tempfile import NamedTemporaryFile
from augur.filter_support.exceptions import FilterException

from augur.dates import (
    ASSERT_ONLY_LESS_SIGNIFICANT_AMBIGUITY_ERROR,
    InvalidDateFormat,
    get_year,
    get_month,
    get_day,
    try_get_numeric_date_min,
    try_get_numeric_date_max,
    get_date_errors,
)
from augur.errors import AugurError
from augur.io.db.sqlite import TabularFileLoaderSQLite, cleanup, ROW_ORDER_COLUMN, sanitize_identifier, chunked_query_to_csv
from augur.io.print import print_err
from augur.utils import read_strains
from augur.filter_support.db.base import DUMMY_COL, FilterBase, FilterCallableReturn, FilterOption
from augur.filter_support.subsample import get_sizes_per_group
from augur.filter_support.output import filter_kwargs_to_str

# internal database globals
# table names
METADATA_TABLE_NAME = 'metadata'
SEQUENCE_INDEX_TABLE_NAME = '_augur_filter_sequence_index'
PRIORITIES_TABLE_NAME = '_augur_filter_priorities'
DATE_TABLE_NAME = '_augur_filter_metadata_date_expanded'
METADATA_FILTER_REASON_TABLE_NAME = '_augur_filter_metadata_filtered_reason'
EXTENDED_FILTERED_TABLE_NAME = '_augur_filter_metadata_filtered_extended'
GROUP_SIZES_TABLE_NAME = '_augur_filter_group_sizes'
OUTPUT_METADATA_TABLE_NAME = '_augur_filter_metadata_output'
# column names
DATE_YEAR_COL = 'year'
DATE_MONTH_COL = 'month'
DATE_DAY_COL = 'day'
NUMERIC_DATE_MIN_COL = 'date_min'
NUMERIC_DATE_MAX_COL = 'date_max'
DATE_ERRORS_COL = 'date_errors'
FILTER_REASON_COL = 'filter'
FILTER_REASON_KWARGS_COL = 'kwargs'
EXCLUDE_COL = 'exclude'
INCLUDE_COL = 'force_include'
GROUP_SIZE_COL = 'size'
PRIORITY_COL = 'priority'
# value for FILTER_REASON_COL with separate logic
SUBSAMPLE_FILTER_REASON = 'subsampling'


class FilterSQLite(FilterBase):
    def __init__(self, args:argparse.Namespace, db_file:str='', in_memory_db=False):
        super().__init__(args)
        self.using_in_memory_db = in_memory_db
        if self.using_in_memory_db:
            # use a singleton connection to an in-memory database https://www.sqlite.org/inmemorydb.html
            self.connection = sqlite3.connect(":memory:")
        else:
            if not db_file:
                tmp_file = NamedTemporaryFile(delete=False)
                db_file = tmp_file.name
            self.db_file = db_file

    def get_db_context(self, **connect_kwargs) -> sqlite3.Connection:
        """Returns a connection to the SQLite database.

        If using an in-memory database, the same connection is returned for all calls to this method.
        Otherwise, a new connection is returned every time."""
        if self.using_in_memory_db:
            return self.connection
        return sqlite3.connect(self.db_file, **connect_kwargs)  # , isolation_level=None

    def db_set_sanitized_identifiers(self):
        """Sets sanitized names for externally sourced identifiers."""
        self.sanitized_metadata_id_column = sanitize_identifier(self.metadata_id_column)
        self.sanitized_date_column = sanitize_identifier(self.date_column)

    def db_create_strain_index(self, table_name:str):
        """Creates a unique index on `metadata_id_column` in a table."""
        index_name = f'idx_{table_name}_id_col'
        with self.get_db_context() as con:
            con.execute(f"""
                CREATE UNIQUE INDEX {index_name}
                ON {table_name} ({self.sanitized_metadata_id_column})
            """)

    def db_load_metadata(self):
        """Loads a metadata file into the database."""
        TabularFileLoaderSQLite(self.args.metadata, self.get_db_context(), METADATA_TABLE_NAME).load()
        try:
            self.db_create_strain_index(METADATA_TABLE_NAME)
        except sqlite3.IntegrityError:
            raise AugurError(f"Duplicate found in '{self.args.metadata}'.")

    def db_load_sequence_index(self):
        """Loads a sequence index file into the database."""
        TabularFileLoaderSQLite(self.args.sequence_index, self.get_db_context(), SEQUENCE_INDEX_TABLE_NAME).load()
        self.db_create_strain_index(SEQUENCE_INDEX_TABLE_NAME)

    def db_get_sequence_index_strains(self):
        """Returns the set of all strains in the sequence index."""
        with self.get_db_context() as con:
            cur = con.execute(f"""
                SELECT {self.sanitized_metadata_id_column}
                FROM {SEQUENCE_INDEX_TABLE_NAME}
            """)
            return {row[0] for row in cur.fetchall()}

    def db_has_date_col(self):
        """Returns a boolean indicating whether `self.date_column` is in the metadata."""
        with self.get_db_context() as con:
            columns = {i[1] for i in con.execute(f'PRAGMA table_info({METADATA_TABLE_NAME})')}
        return (self.date_column in columns)

    def db_create_date_table(self):
        """Creates an intermediate date table from the metadata table.

        Contains the strain column, original date column, and these extracted columns:
        - `DATE_YEAR_COL`: Extracted year (int or `NULL`)
        - `DATE_MONTH_COL`: Extracted month (int or `NULL`)
        - `DATE_DAY_COL`: Extracted day (int or `NULL`)
        - `DATE_MIN_COL`: Exact date, minimum if ambiguous
        - `DATE_MAX_COL`: Exact date, maximum if ambiguous
        """
        with self.get_db_context() as con:
            if self.has_date_col:
                # TODO: handle numeric dates for year/month/day
                con.create_function(get_year.__name__, 1, get_year)
                con.create_function(get_month.__name__, 1, get_month)
                con.create_function(get_day.__name__, 1, get_day)
                con.create_function(try_get_numeric_date_min.__name__, 1, try_get_numeric_date_min)
                con.create_function(try_get_numeric_date_max.__name__, 1, try_get_numeric_date_max)
                con.create_function(get_date_errors.__name__, 1, get_date_errors)
                con.execute(f"""CREATE TABLE {DATE_TABLE_NAME} AS
                    SELECT
                        {self.sanitized_metadata_id_column},
                        {self.sanitized_date_column},
                        {get_year.__name__}({self.sanitized_date_column}) AS {DATE_YEAR_COL},
                        {get_month.__name__}({self.sanitized_date_column}) AS {DATE_MONTH_COL},
                        {get_day.__name__}({self.sanitized_date_column}) AS {DATE_DAY_COL},
                        {try_get_numeric_date_min.__name__}({self.sanitized_date_column}) AS {NUMERIC_DATE_MIN_COL},
                        {try_get_numeric_date_max.__name__}({self.sanitized_date_column}) AS {NUMERIC_DATE_MAX_COL},
                        {get_date_errors.__name__}({self.sanitized_date_column}) AS {DATE_ERRORS_COL}
                    FROM {METADATA_TABLE_NAME}
                """)
                self._validate_date_table()
            else:
                # create placeholder table for later JOINs
                con.execute(f"""CREATE TABLE {DATE_TABLE_NAME} AS
                    SELECT
                        {self.sanitized_metadata_id_column},
                        '' AS {DATE_YEAR_COL},
                        '' AS {DATE_MONTH_COL},
                        '' AS {DATE_DAY_COL},
                        '' AS {NUMERIC_DATE_MIN_COL},
                        '' AS {NUMERIC_DATE_MAX_COL}
                    FROM {METADATA_TABLE_NAME}
                """)
        self.db_create_strain_index(DATE_TABLE_NAME)

    def _validate_date_table(self):
        """Validate dates in `DATE_TABLE_NAME`.

        Internally runs a query for invalid dates, i.e. rows where:
        1. date was specified (not null or empty string)
        2. date failed assert_only_less_significant_ambiguity check

        Raises
        ------
        :class:`InvalidDateFormat`
        """
        max_results = 3 # limit length of error message
        with self.get_db_context() as con:
            cur = con.execute(f"""
                SELECT CAST({self.sanitized_date_column} AS TEXT)
                FROM {DATE_TABLE_NAME}
                WHERE NOT ({self.sanitized_date_column} IS NULL OR {self.sanitized_date_column} = '')
                    AND ({DATE_ERRORS_COL} = '{ASSERT_ONLY_LESS_SIGNIFICANT_AMBIGUITY_ERROR}')
                LIMIT {max_results}
            """)
            invalid_dates = [repr(row[0]) for row in cur.fetchall()]
        if invalid_dates:
            raise InvalidDateFormat(f"Some dates have an invalid format (showing at most {max_results}): {','.join(invalid_dates)}.\n"
                + "If year contains ambiguity, month and day must also be ambiguous.\n"
                + "If month contains ambiguity, day must also be ambiguous.")

    def filter_by_exclude_all(self) -> FilterCallableReturn:
        """Exclude all strains regardless of the given metadata content.

        This is a placeholder function that can be called as part of a generalized
        loop through all possible functions.

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        expression = 'True'
        parameters = {}
        return expression, parameters

    def filter_by_exclude_strains(self, exclude_file) -> FilterCallableReturn:
        """Exclude the given set of strains from the given metadata.

        Parameters
        ----------
        exclude_file : str
            Filename with strain names to exclude from the given metadata

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        excluded_strains = read_strains(exclude_file)
        excluded_strains = [f"'{strain}'" for strain in excluded_strains]
        expression = f"""
            {self.sanitized_metadata_id_column} IN ({','.join(excluded_strains)})
        """
        parameters = {}
        return expression, parameters

    def parse_filter_query(self, query):
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
        callable :
            Operator function to test equality or non-equality of values
        str :
            Value of column to query

        >>> FilterSQLite.parse_filter_query(FilterSQLite, "property=value")
        ('property', '=', 'value')
        >>> FilterSQLite.parse_filter_query(FilterSQLite, "property!=value")
        ('property', '!=', 'value')

        """
        column, value = re.split(r'!?=', query)
        op = '='
        if "!=" in query:
            op = '!='

        return column, op, value

    def filter_by_exclude_where(self, exclude_where) -> FilterCallableReturn:
        """Exclude all strains from the given metadata that match the given exclusion query.

        Unlike pandas query syntax, exclusion queries should follow the pattern of
        `"property=value"` or `"property!=value"`. Additionally, this filter treats
        all values like lowercase strings, so we convert all values to strings first
        and then lowercase them before testing the given query.

        Parameters
        ----------
        exclude_where : str
            Filter query used to exclude strains

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        column, op, value = self.parse_filter_query(exclude_where)
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {METADATA_TABLE_NAME}
                WHERE {METADATA_TABLE_NAME}.{sanitize_identifier(column)} {op} :value
            )
        """
        parameters = {'value': value}
        return expression, parameters

    def filter_by_query(self, query) -> FilterCallableReturn:
        """Filter by any valid SQL expression on the metadata.

        Strains that do *not* match the query will be excluded.

        Parameters
        ----------
        query : str
            SQL expression used to exclude strains

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        # NOT query to exclude all that do not match
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {METADATA_TABLE_NAME}
                WHERE NOT ({query})
            )
        """
        parameters = {}
        return expression, parameters

    def filter_by_ambiguous_date(self, ambiguity="any") -> FilterCallableReturn:
        """Filter metadata in the given pandas DataFrame where values in the given date
        column have a given level of ambiguity.

        Determine ambiguity hierarchically such that, for example, an ambiguous
        month implicates an ambiguous day even when day information is available.

        Parameters
        ----------
        ambiguity : str
            Level of date ambiguity to filter metadata by

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        if ambiguity == 'year':
            expression = f"""
                {self.sanitized_metadata_id_column} IN (
                    SELECT {self.sanitized_metadata_id_column}
                    FROM {DATE_TABLE_NAME}
                    WHERE {DATE_YEAR_COL} IS NULL
                )
            """
        elif ambiguity == 'month':
            expression = f"""
                {self.sanitized_metadata_id_column} IN (
                    SELECT {self.sanitized_metadata_id_column}
                    FROM {DATE_TABLE_NAME}
                    WHERE {DATE_MONTH_COL} IS NULL OR {DATE_YEAR_COL} IS NULL
                )
            """
        elif ambiguity == 'day' or ambiguity == 'any':
            expression = f"""
                {self.sanitized_metadata_id_column} IN (
                    SELECT {self.sanitized_metadata_id_column}
                    FROM {DATE_TABLE_NAME}
                    WHERE {DATE_DAY_COL} IS NULL OR {DATE_MONTH_COL} IS NULL OR {DATE_YEAR_COL} IS NULL
                )
            """
        parameters = {}
        return expression, parameters

    def filter_by_min_date(self, min_date:float) -> FilterCallableReturn:
        """Filter metadata by minimum date.

        Parameters
        ----------
        min_date : float
            Minimum date

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {DATE_TABLE_NAME}
                WHERE {NUMERIC_DATE_MAX_COL} < :min_date OR {NUMERIC_DATE_MIN_COL} IS NULL
            )
        """
        parameters = {'min_date': min_date}
        return expression, parameters

    def filter_by_max_date(self, max_date:float) -> FilterCallableReturn:
        """Filter metadata by maximum date.

        Parameters
        ----------
        max_date : float
            Maximum date

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {DATE_TABLE_NAME}
                WHERE {NUMERIC_DATE_MIN_COL} > :max_date OR {NUMERIC_DATE_MAX_COL} IS NULL
            )
        """
        parameters = {'max_date': max_date}
        return expression, parameters

    def filter_by_sequence_index(self) -> FilterCallableReturn:
        """Filter metadata by presence of corresponding entries in a given sequence
        index. This filter effectively intersects the strain ids in the metadata and
        sequence index.

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        expression = f"""
            {self.sanitized_metadata_id_column} NOT IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {SEQUENCE_INDEX_TABLE_NAME}
            )
        """
        parameters = {}
        return expression, parameters

    def filter_by_sequence_length(self, min_length=0) -> FilterCallableReturn:
        """Filter metadata by sequence length from a given sequence index.

        Parameters
        ----------
        min_length : int
            Minimum number of standard nucleotide characters (A, C, G, or T) in each sequence

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {SEQUENCE_INDEX_TABLE_NAME}
                WHERE A+C+G+T < :min_length
            )
        """
        parameters = {'min_length': min_length}
        return expression, parameters

    def filter_by_non_nucleotide(self) -> FilterCallableReturn:
        """Filter metadata for strains with invalid nucleotide content.

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {SEQUENCE_INDEX_TABLE_NAME}
                WHERE invalid_nucleotides != 0
            )
        """
        parameters = {}
        return expression, parameters

    def force_include_strains(self, include_file) -> FilterCallableReturn:
        """Include strains in the given text file from the given metadata.

        Parameters
        ----------
        include_file : str
            Filename with strain names to include from the given metadata

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        included_strains = read_strains(include_file)
        included_strains = [f"'{strain}'" for strain in included_strains]
        expression = f"""
            {self.sanitized_metadata_id_column} IN ({','.join(included_strains)})
        """
        parameters = {}
        return expression, parameters

    def force_include_where(self, include_where) -> FilterCallableReturn:
        """Include all strains from the given metadata that match the given query.

        Unlike pandas query syntax, inclusion queries should follow the pattern of
        `"property=value"` or `"property!=value"`. Additionally, this filter treats
        all values like lowercase strings, so we convert all values to strings first
        and then lowercase them before testing the given query.

        Parameters
        ----------
        include_where : str
            Filter query used to include strains

        Returns
        -------
        (str, dict):
            str: expression for SQL query `WHERE` clause
            dict: named parameters used in the expression, if any
        """
        column, op, value = self.parse_filter_query(include_where)
        expression = f"""
            {self.sanitized_metadata_id_column} IN (
                SELECT {self.sanitized_metadata_id_column}
                FROM {METADATA_TABLE_NAME}
                WHERE {METADATA_TABLE_NAME}.{sanitize_identifier(column)} {op} :value
            )
        """
        parameters = {'value': value}
        return expression, parameters

    def db_create_filter_reason_table(self, exclude_by:List[FilterOption], include_by:List[FilterOption]):
        """Creates an intermediate table for filter reason.
        
        Applies exclusion and force-inclusion rules to filter strains from the metadata.

        Parameters
        ----------
        exclude_by : list
            A list of filter expressions for SQL query `WHERE` clause
        include_by : list
            A list of filter expressions for SQL query `WHERE` clause
        """
        with self.get_db_context() as con:
            con.execute(f"""CREATE TABLE {METADATA_FILTER_REASON_TABLE_NAME} AS
                SELECT
                    {self.sanitized_metadata_id_column},
                    FALSE AS {EXCLUDE_COL},
                    FALSE AS {INCLUDE_COL},
                    NULL AS {FILTER_REASON_COL},
                    NULL AS {FILTER_REASON_KWARGS_COL}
                FROM {METADATA_TABLE_NAME}
            """)
        self.db_create_strain_index(METADATA_FILTER_REASON_TABLE_NAME)
        # note: consider JOIN vs subquery if performance issues https://stackoverflow.com/q/3856164
        self.db_apply_exclusions(exclude_by)
        self.db_apply_force_inclusions(include_by)

    def db_apply_exclusions(self, exclude_by:List[FilterOption]):
        """Updates the filter reason table with exclusion rules."""
        for exclude_function, kwargs in exclude_by:
            where_expression, where_parameters = exclude_function(**kwargs)
            sql = f"""
                UPDATE {METADATA_FILTER_REASON_TABLE_NAME}
                SET
                    {EXCLUDE_COL} = TRUE,
                    {FILTER_REASON_COL} = :filter_reason,
                    {FILTER_REASON_KWARGS_COL} = :filter_reason_kwargs
                WHERE {where_expression}
            """
            sql_parameters = {
                'filter_reason': exclude_function.__name__,
                'filter_reason_kwargs': filter_kwargs_to_str(kwargs)
            }
            sql_parameters = {**sql_parameters, **where_parameters}
            try:
                with self.get_db_context() as con:
                    con.execute(sql, sql_parameters)
            except sqlite3.OperationalError as sql_e:
                if str(sql_e).startswith('no such column'):
                    raise FilterException(sql_e) from sql_e

    def db_apply_force_inclusions(self, include_by:List[FilterOption]):
        """Updates the filter reason table with force-inclusion rules."""
        for include_function, kwargs in include_by:
            where_expression, where_parameters = include_function(**kwargs)
            sql = f"""
                UPDATE {METADATA_FILTER_REASON_TABLE_NAME}
                SET
                    {INCLUDE_COL} = TRUE,
                    {FILTER_REASON_COL} = :filter_reason,
                    {FILTER_REASON_KWARGS_COL} = :filter_reason_kwargs
                WHERE {where_expression}
            """
            sql_parameters = {
                'filter_reason': include_function.__name__,
                'filter_reason_kwargs': filter_kwargs_to_str(kwargs)
            }
            sql_parameters = {**sql_parameters, **where_parameters}
            try:
                with self.get_db_context() as con:
                    con.execute(sql, sql_parameters)
            except sqlite3.OperationalError as sql_e:
                if str(sql_e).startswith('no such column'):
                    raise FilterException(sql_e) from sql_e

    def db_create_output_table(self):
        """Creates a final intermediate table to be used for output.

        The table is a subset of the original metadata table containing rows that pass all:
        1. exclusion rules
        2. force-inclusion rules
        3. subsampling
        """
        with self.get_db_context() as con:
            con.execute(f"""CREATE TABLE {OUTPUT_METADATA_TABLE_NAME} AS
                SELECT m.* FROM {METADATA_TABLE_NAME} m
                JOIN {METADATA_FILTER_REASON_TABLE_NAME} f
                    USING ({self.sanitized_metadata_id_column})
                WHERE NOT f.{EXCLUDE_COL} OR f.{INCLUDE_COL}
            """)

    def db_get_counts_per_group(self, group_by_cols:List[str]) -> List[int]:
        """
        Returns
        -------
        list[int]
            List of counts per group.
        """
        sanitized_group_by_cols = [sanitize_identifier(col) for col in group_by_cols]
        with self.get_db_context() as con:
            cur = con.execute(f"""
                SELECT {','.join(sanitized_group_by_cols)}, COUNT(*)
                FROM {EXTENDED_FILTERED_TABLE_NAME}
                GROUP BY {','.join(sanitized_group_by_cols)}
            """)
            return [row[-1] for row in cur.fetchall()]

    def db_get_filtered_strains_count(self) -> int:
        """Returns the number of metadata strains that pass all filter rules.

        Note: this can return a different number before and after subsampling.
        """
        with self.get_db_context() as con:
            return con.execute(f"""
                SELECT COUNT(*)
                FROM {METADATA_FILTER_REASON_TABLE_NAME}
                WHERE NOT {EXCLUDE_COL} OR {INCLUDE_COL}
            """).fetchone()[0]

    def db_update_filter_reason_table_with_subsampling(self, group_by_cols:List[str]):
        """Subsamples filtered metadata and updates the filter reason table."""
        sanitized_group_by_cols = [sanitize_identifier(col) for col in group_by_cols]
        # create a SQL query for strains to subsample for
        where_conditions = [f'group_i <= {GROUP_SIZE_COL}']
        for col in sanitized_group_by_cols:
            where_conditions.append(f'{col} IS NOT NULL')
        # `ORDER BY ... NULLS LAST` is unsupported for SQLite <3.30.0 so `CASE ... IS NULL` is a workaround
        # ref https://stackoverflow.com/a/12503284
        query_for_subsampled_strains = f"""
            SELECT {self.sanitized_metadata_id_column}
            FROM (
                SELECT {self.sanitized_metadata_id_column}, {','.join(sanitized_group_by_cols)}, ROW_NUMBER() OVER (
                    PARTITION BY {','.join(sanitized_group_by_cols)}
                    ORDER BY (CASE WHEN {PRIORITY_COL} IS NULL THEN 1 ELSE 0 END), {PRIORITY_COL} DESC
                ) AS group_i
                FROM {EXTENDED_FILTERED_TABLE_NAME}
            )
            JOIN {GROUP_SIZES_TABLE_NAME} USING({','.join(sanitized_group_by_cols)})
            WHERE {' AND '.join(where_conditions)}
        """
        # update filter reason table
        with self.get_db_context() as con:
            con.execute(f"""
                UPDATE {METADATA_FILTER_REASON_TABLE_NAME}
                SET
                    {EXCLUDE_COL} = TRUE,
                    {FILTER_REASON_COL} = '{SUBSAMPLE_FILTER_REASON}'
                WHERE NOT {EXCLUDE_COL} AND {self.sanitized_metadata_id_column} NOT IN (
                    {query_for_subsampled_strains}
                )
            """)

    def db_load_priorities_table(self):
        """Loads a priorities file into the database."""
        try:
            TabularFileLoaderSQLite(self.args.priority, self.get_db_context(), PRIORITIES_TABLE_NAME,
                    header=False, names=[self.metadata_id_column, PRIORITY_COL]).load()
        except ValueError as e:
            raise ValueError(f"Failed to parse priority file {self.args.priority}.") from e
        self.db_create_strain_index(PRIORITIES_TABLE_NAME)

    def db_generate_priorities_table(self, seed:int=None):
        """Generates a priorities table with random priorities."""
        # Use pandas/numpy since random seeding is not possible with SQLite:
        # https://stackoverflow.com/a/24394275
        df_priority = pd.read_sql(f"""
                SELECT {self.sanitized_metadata_id_column}
                FROM {METADATA_FILTER_REASON_TABLE_NAME}
                WHERE NOT {EXCLUDE_COL} OR {INCLUDE_COL}
            """, self.get_db_context())
        rng = np.random.default_rng(seed)
        df_priority[PRIORITY_COL] = rng.random(len(df_priority))
        df_priority.to_sql(PRIORITIES_TABLE_NAME, self.get_db_context(), index=False)
        self.db_create_strain_index(PRIORITIES_TABLE_NAME)

    def db_get_metadata_cols(self) -> Set[str]:
        """Returns a set of metadata column names."""
        with self.get_db_context() as con:
            return {i[1] for i in con.execute(f'PRAGMA table_info({METADATA_TABLE_NAME})')}

    def db_create_extended_filtered_metadata_table(self, group_by_cols:List[str]):
        """Creates a new table with rows as filtered metadata, with the following columns:

        1. Metadata ID column
        2. Group-by columns
        2. Extracted date columns (`DATE_YEAR_COL`, `DATE_MONTH_COL`)
        3. `PRIORITY_COL`
        4. `dummy` containing the same value in all rows, used when no group-by columns are provided
        """
        if DATE_YEAR_COL in group_by_cols and DATE_YEAR_COL in self.db_get_metadata_cols():
            print_err(f"WARNING: `--group-by year` uses the generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.")
        if DATE_MONTH_COL in group_by_cols and DATE_MONTH_COL in self.db_get_metadata_cols():
            print_err(f"WARNING: `--group-by month` uses the generated month value from the 'date' column. The custom 'month' column in the metadata is ignored for grouping purposes.")

        # Start with an empty string in case this isn't needed (e.g. when group_by_cols=['year', 'month']).
        sanitized_group_by_cols_for_select = ''
        if group_by_cols:
            # Copy to avoid modifying the original list.
            group_by_cols_copy = list(group_by_cols)

            # Ignore extracted date columns - those are added directly from the date table.
            if DATE_YEAR_COL in group_by_cols:
                group_by_cols_copy.remove(DATE_YEAR_COL)
            if DATE_MONTH_COL in group_by_cols:
                group_by_cols_copy.remove(DATE_MONTH_COL)

            if group_by_cols_copy:
                # Prefix columns with the metadata table alias defined in the SQL query further down.
                sanitized_group_by_cols_for_select = ','.join(f'm.{sanitize_identifier(col)}' for col in group_by_cols_copy)
                # Add an extra comma for valid SQL.
                sanitized_group_by_cols_for_select += ','

        # LEFT OUTER JOIN the priorities table since a default inner join would
        # drop strains without a priority.
        with self.get_db_context() as con:
            con.execute(f"""CREATE TABLE {EXTENDED_FILTERED_TABLE_NAME} AS
                SELECT
                    m.{self.sanitized_metadata_id_column},
                    {sanitized_group_by_cols_for_select}
                    d.{DATE_YEAR_COL}, d.{DATE_MONTH_COL},
                    p.{PRIORITY_COL},
                    TRUE AS {DUMMY_COL}
                FROM {METADATA_TABLE_NAME} m
                JOIN {METADATA_FILTER_REASON_TABLE_NAME} f ON (m.{self.sanitized_metadata_id_column} = f.{self.sanitized_metadata_id_column})
                    AND (NOT f.{EXCLUDE_COL} OR f.{INCLUDE_COL})
                JOIN {DATE_TABLE_NAME} d USING ({self.sanitized_metadata_id_column})
                LEFT OUTER JOIN {PRIORITIES_TABLE_NAME} p USING ({self.sanitized_metadata_id_column})
            """)

    def db_create_group_sizes_table(self, group_by_cols:List[str], sequences_per_group:float):
        """Creates an intermediate table for group sizes."""
        sanitized_group_by_cols = [sanitize_identifier(col) for col in group_by_cols]

        # Get a DataFrame with one row per unique group.
        df_groups = pd.read_sql_query(f"""
                SELECT {','.join(sanitized_group_by_cols)}
                FROM {EXTENDED_FILTERED_TABLE_NAME}
                GROUP BY {','.join(sanitized_group_by_cols)}
            """, self.get_db_context())

        # Extend the above DataFrame with a new column that indicates the number of samples to include in each group.
        df_sizes = get_sizes_per_group(df_groups, GROUP_SIZE_COL, sequences_per_group, random_seed=self.args.subsample_seed)

        # Write back to a SQL table.
        df_sizes.to_sql(GROUP_SIZES_TABLE_NAME, self.get_db_context())

    def db_output_strains(self):
        """Writes the final output strains to a file, one strain per line.
        Row order from the original metadata file is retained.
        """
        query = f"""
            SELECT {self.sanitized_metadata_id_column}
            FROM {OUTPUT_METADATA_TABLE_NAME}
            ORDER BY {ROW_ORDER_COLUMN}
        """
        chunked_query_to_csv(
            con=self.get_db_context(),
            query=query,
            path=self.args.output_strains,
            chunksize=10000,
            header=False,
            index=False,
        )

    def db_output_metadata(self):
        """Writes the final output strains to a file, one strain per line,
        along with the original metadata.
        Row order from the original metadata file is retained.
        """
        query = f"""
            SELECT *
            FROM {OUTPUT_METADATA_TABLE_NAME}
            ORDER BY {ROW_ORDER_COLUMN}
        """
        chunked_query_to_csv(
            con=self.get_db_context(),
            query=query,
            path=self.args.output_metadata,
            chunksize=10000,
            header=True,
            index=False,
            columns_to_exclude=[ROW_ORDER_COLUMN],
            sep='\t',
        )

    def db_output_log(self):
        """Writes a file explaining the reason for excluded or force-included strains.

        This file has the following columns:
        1. Original strain column
        2. Name of the filter function responsible for inclusion/exclusion
        3. Arguments given to the filter function
        """
        query = f"""
            SELECT {self.sanitized_metadata_id_column}, {FILTER_REASON_COL}, {FILTER_REASON_KWARGS_COL}
            FROM {METADATA_FILTER_REASON_TABLE_NAME}
            WHERE {FILTER_REASON_COL} IS NOT NULL
        """
        chunked_query_to_csv(
            con=self.get_db_context(),
            query=query,
            path=self.args.output_log,
            chunksize=10000,
            header=True,
            index=False,
            sep='\t',
        )

    def db_get_metadata_strains(self) -> Set[str]:
        """Returns a set of all strain names from the original metadata."""
        with self.get_db_context() as con:
            cur = con.execute(f"""
                SELECT {self.sanitized_metadata_id_column}
                FROM {METADATA_TABLE_NAME}
            """)
            return {row[0] for row in cur.fetchall()}

    def db_get_strains_passed(self) -> Set[str]:
        """Returns a set of all strain names that passed filtering and subsampling."""
        with self.get_db_context() as con:
            cur = con.execute(f"""
                SELECT {self.sanitized_metadata_id_column}
                FROM {METADATA_FILTER_REASON_TABLE_NAME}
                WHERE NOT {EXCLUDE_COL} OR {INCLUDE_COL}
            """)
            return {row[0] for row in cur.fetchall()}

    def db_get_num_metadata_strains(self) -> int:
        """Returns the number of strains in the original metadata."""
        with self.get_db_context() as con:
            return con.execute(f"""
                SELECT COUNT(*)
                FROM {METADATA_TABLE_NAME}
            """).fetchone()[0]

    def db_get_num_excluded_by_subsampling(self) -> int:
        """Returns the number of strains excluded by subsampling."""
        with self.get_db_context() as con:
            return con.execute(f"""
                SELECT COUNT(*)
                FROM {METADATA_FILTER_REASON_TABLE_NAME}
                WHERE {FILTER_REASON_COL} = '{SUBSAMPLE_FILTER_REASON}'
            """).fetchone()[0]

    def db_get_filter_counts(self) -> List[Tuple[str, str, int]]:
        """Returns a list of tuples for each filter function that had an effect.

        Each tuple has:
        1. Name of the filter function
        2. Arguments given to the filter function
        3. Number of strains included/excluded by the filter function
        """
        with self.get_db_context() as con:
            return con.execute(f"""
                SELECT {FILTER_REASON_COL}, {FILTER_REASON_KWARGS_COL}, COUNT(*)
                FROM {METADATA_FILTER_REASON_TABLE_NAME}
                WHERE {FILTER_REASON_COL} IS NOT NULL
                    AND {FILTER_REASON_COL} != '{SUBSAMPLE_FILTER_REASON}'
                GROUP BY {FILTER_REASON_COL}, {FILTER_REASON_KWARGS_COL}
            """).fetchall()

    def db_cleanup(self):
        """Deletes any temporary files."""
        if not self.using_in_memory_db:
            cleanup(self.db_file)
