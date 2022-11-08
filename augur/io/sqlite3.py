import os
import pandas as pd
import random
import re
import string
import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from typing import Any, Dict, Iterable, List

from augur.errors import AugurError


class Sqlite3DatabaseError(AugurError):
    pass


class Sqlite3Database:
    """Represents a SQLite3 database.

    Provides useful methods to manage connections and general I/O.
    """

    def __init__(self, file: str):

        self.file = file
        """Database file."""

        self.connection: sqlite3.Connection = None
        """SQLite3 database connection."""

    def connect(self, **connect_kwargs) -> sqlite3.Connection:
        """Return a new connection to the SQLite database."""
        return sqlite3.connect(self.file, **connect_kwargs)

    def __enter__(self):
        self.connection = self.connect().__enter__()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.connection.__exit__(exc_type, exc_value, exc_traceback)

    def tables(self):
        """Yield table names in the database."""
        for row in self.execute(f"SELECT name FROM sqlite_master where type='table'"):
            yield row[0]

    def columns(self, table_name: str):
        """Yield columns names from a table."""
        for row in self.execute(f"SELECT name FROM pragma_table_info({sanitize_identifier(table_name)})"):
            yield row[0]

    def execute(self, *args):
        """Execute under a connection context."""
        if self.connection:
            return self.connection.execute(*args)
        else:
            with self.connect() as con:
                return con.execute(*args)

    def executemany(self, *args):
        """Executemany under a connection context."""
        if self.connection:
            return self.connection.executemany(*args)
        else:
            with self.connect() as con:
                return con.executemany(*args)

    def check_table_existence(self, table_name: str):
        """Raise an error if the table does not exist."""
        if table_name not in self.tables():
            raise Sqlite3DatabaseError(f"Table '{table_name}' does not exist.")

    def create_table(self, table_name: str, columns: List[str]):
        """Create a table with all columns as TEXT type."""
        try:
            self.execute(f"""
                CREATE TABLE {sanitize_identifier(table_name)} (
                    {','.join(f'{sanitize_identifier(column)} TEXT' for column in columns)}
                )
            """)
        except sqlite3.OperationalError as e:
            raise Sqlite3DatabaseError(f'Failed to create table: {e}')

    def insert(self, table_name: str, column_names: List[str], rows: Iterable[Dict[str, Any]]):
        """Insert rows into a table."""
        # This intermediate OrderedDict serves two purposes:
        # 1. Generates valid placeholder names. Column names cannot be used directly
        #    since placeholder names cannot be quote-sanitized in the insert
        #    statement, and column names are user-defined.
        # 2. Ensures the list of column names and the list of placeholders in the insert statement are 1:1.
        column_placeholder_mapping = OrderedDict({col: _generate_placeholder(col) for col in column_names})

        insert_statement = f"""
            INSERT INTO {sanitize_identifier(table_name)}
            ({','.join([sanitize_identifier(col) for col in column_placeholder_mapping.keys()])})
            VALUES ({','.join([f':{placeholder}' for placeholder in column_placeholder_mapping.values()])})
        """
        rows_with_placeholders = (
            {column_placeholder_mapping[col]: value
                for col, value in row.items()
                if col
            }
            for row in rows
        )
        try:
            self.executemany(insert_statement, rows_with_placeholders)
        except sqlite3.ProgrammingError as e:
            raise Sqlite3DatabaseError("Failed to import rows.") from e

    def table_to_csv(self, table_name: str, path: str, header: bool = True, chunksize: int = 10000, **to_csv_kwargs):
        """Output a table's contents to a tabular file.

        Parameters
        ----------
        table_name
            Table to output.
        path
            Path to the output file.
        header
            Write out the column names.
        chunksize
            Maximum number of rows to load into memory at once.
        to_csv_kwargs
            Additional keyword arguments passed to pandas.to_csv().
        """
        self.check_table_existence(table_name)
        query = f"""
            SELECT * FROM {sanitize_identifier(table_name)}
        """

        with self.connect() as con:
            # Get results in chunks so they are not loaded into memory all at once.
            df_chunks = pd.read_sql_query(query, con, chunksize=chunksize)

            # Overwrite the file if it exists.
            mode = 'w'

            for df_chunk in df_chunks:
                df_chunk.to_csv(path, mode=mode, header=header, **to_csv_kwargs)

                # For subsequent chunks, use append mode and don't include the header.
                mode = 'a'
                header = False

    def table_to_fasta(self, table_name: str, path: str, id_column: str, sequence_column: str, chunksize: int = 10000):
        """Output two columns from a table as a FASTA file.
        
        Parameters
        ----------
        table_name
            Table to output.
        path
            Path to the output file.
        id_column
            Column to use for SeqID header lines.
        sequence_column
            Column to use as the actual sequence following header lines.
        chunksize
            Maximum number of rows to load into memory at once.
        """
        self.check_table_existence(table_name)
        columns = set(self.columns(table_name))
        if id_column not in columns:
            raise Sqlite3DatabaseError(f"Sequence ID column '{id_column}' missing from table '{table_name}'.")
        if sequence_column not in columns:
            raise Sqlite3DatabaseError(f"Sequence column '{sequence_column}' missing from table '{table_name}'.")

        query = f"""
            SELECT {sanitize_identifier(id_column)},
                {sanitize_identifier(sequence_column)}
            FROM {sanitize_identifier(table_name)}
        """

        with self.connect() as con, open(path, "w") as out_file:
            # Get results in chunks so they are not loaded into memory all at once.
            df_chunks = pd.read_sql_query(query, con, chunksize=chunksize)

            for df_chunk in df_chunks:
                sequences = df_chunk.apply(lambda row: SeqRecord(id=row[id_column], seq=Seq(row[sequence_column]), description=""), axis=1)
                if not sequences.empty:
                    SeqIO.write(sequences, out_file, "fasta")

    def delete(self):
        """Remove the database file if present."""
        try:
            os.remove(self.file)
        except FileNotFoundError:
            pass


def sanitize_identifier(identifier: str):
    """Sanitize a SQLite identifier."""
    # Note: we can (and probably should) do more here.
    # https://stackoverflow.com/a/6701665

    # Escape existing double quotes
    identifier = identifier.replace('"', '""')

    # Wrap inside double quotes
    return f'"{identifier}"'


def _generate_placeholder(identifier: str):
    """Generate a placeholder name that doesn't need escaping.

    More on placeholders: https://www.sqlite.org/lang_expr.html#varparam
    """
    # Remove everything but alphanumeric characters and underscores
    stripped_identifier = re.sub(r'\W+', '', identifier)

    # Generate a random suffix to avoid collisions
    random_suffix = ''.join(random.choices(string.ascii_letters + string.digits, k=10))

    return f'{stripped_identifier}_{random_suffix}'
