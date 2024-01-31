import csv
import random
import re
import string
import sqlite3
from contextlib import AbstractContextManager
from typing import Any, Dict, Iterable
from urllib.parse import urlencode
from .file import open_file


class Sqlite3Database(AbstractContextManager):
    """Represents a SQLite3 database.

    This class should be used as a context manager where the runtime context of
    an instance reflects a connection to the database.
    """

    def __init__(self, file: str, **connect_uri_params):
        """
        Parameters
        ----------
        file
            Database file
        connect_uri_params
            Parameters passed to sqlite3.connect() in the form of URI parameters.
            Examples: https://docs.python.org/3/library/sqlite3.html#how-to-work-with-sqlite-uris
        """

        self.file = file
        """Database file."""

        self.connection: sqlite3.Connection = None
        """SQLite3 database connection."""

        self.connect_uri_params = connect_uri_params
        """Parameters passed to sqlite3.connect() in the form of URI parameters."""

        # Default to opening the database in read-only mode
        if 'mode' not in self.connect_uri_params:
            self.connect_uri_params['mode'] = 'ro'

    def _connect(self, **connect_uri_params):
        """Return a new connection to the SQLite database.

        Intended for use with the context manager."""
        encoded_params = urlencode(connect_uri_params)
        return sqlite3.connect(f'file:{self.file}?{encoded_params}', uri=True)

    def __enter__(self):
        self.connection = self._connect(**self.connect_uri_params).__enter__()

        # Allow access by column name.
        self.connection.row_factory = sqlite3.Row

        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.connection.__exit__(exc_type, exc_value, exc_traceback)

    def tables(self):
        """Yield table names in the database."""
        for row in self.connection.execute(f"SELECT name FROM sqlite_master WHERE type='table'"):
            yield str(row["name"])

    def columns(self, table: str):
        """Yield columns names from a table."""
        self._assert_table_existence(table)

        for row in self.connection.execute(f"SELECT name FROM pragma_table_info({sanitize_identifier(table)})"):
            yield str(row["name"])

    def _assert_table_existence(self, table: str):
        """Assert that the table exists."""
        assert table in self.tables()

    def create_table(self, table: str, columns: Iterable[str]):
        """Create a table with all columns having the TEXT type affinity."""

        # Create table with TEXT type affinity for all columns.
        # FIXME: STRICT requires SQLite version 3.37.0. Do we actually need it?
        try:
            self.connection.execute(f"""
                CREATE TABLE {sanitize_identifier(table)} (
                    {','.join(f'{sanitize_identifier(column)} TEXT' for column in columns)}
                ) STRICT
            """)
        except sqlite3.OperationalError as e:
            raise Exception(f'Failed to create table.') from e

    def insert(self, table: str, columns: Iterable[str], rows: Iterable[Dict[str, Any]]):
        """Insert rows into a table."""
        self._assert_table_existence(table)

        # This intermediate dict serves two purposes:
        # 1. Generates valid placeholder names. Column names cannot be used
        #    directly since placeholder names cannot be quote-sanitized in the
        #    insert statement, and column names are user-defined.
        # 2. Ensures the list of column names and the list of placeholders in
        #    the insert statement are 1:1.
        column_placeholder_mapping = {column: _generate_placeholder(column) for column in columns}

        insert_statement = f"""
            INSERT INTO {sanitize_identifier(table)}
            ({','.join([sanitize_identifier(column) for column in column_placeholder_mapping.keys()])})
            VALUES ({','.join([f':{placeholder}' for placeholder in column_placeholder_mapping.values()])})
        """
        rows_with_placeholders = (
            {column_placeholder_mapping[key]: value
                for key, value in row.items()
                if key in columns
            }
            for row in rows
        )
        try:
            self.connection.executemany(insert_statement, rows_with_placeholders)
        except sqlite3.ProgrammingError as e:
            raise Exception("Failed to insert rows.") from e

    def query_to_file(self, query: str, path: str, header: bool = True, delimiter: str = '\t'):
        """Query the database and write results to a tabular file.

        Parameters
        ----------
        query
            SQLite query string.
        path
            Path to the output file.
        header
            Write out the column names.
        delimiter
            Delimiter between columns.
        """
        with open_file(path, mode="w") as output_file:
            cursor = self.connection.cursor()
            cursor.execute(query)
            writer = csv.writer(output_file, delimiter=delimiter, lineterminator='\n')
            if header:
                column_names = [column[0] for column in cursor.description]
                writer.writerow(column_names)
            for row in cursor:
                writer.writerow(row)


def sanitize_identifier(identifier: str):
    """Sanitize a SQLite identifier.

    Note: We can (and probably should) do more here¹. However, column names in
    the database are used in the output, which should be as close as possible
    to the input names.
    ¹ https://stackoverflow.com/a/6701665
    """

    # Escape existing double quotes
    identifier = identifier.replace('"', '""')

    # Wrap inside double quotes
    return f'"{identifier}"'


def _generate_placeholder(identifier: str):
    """Generate a placeholder¹ name that doesn't need escaping.

    ¹ https://www.sqlite.org/lang_expr.html#varparam
    """
    # Remove everything but alphanumeric characters and underscores
    stripped_identifier = re.sub(r'\W+', '', identifier)

    # Generate a random suffix to avoid collisions
    random_suffix = ''.join(random.choices(string.ascii_letters + string.digits, k=10))

    return f'{stripped_identifier}_{random_suffix}'
