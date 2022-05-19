import os
import pandas as pd
import sqlite3
from typing import List

from . import TabularFileLoaderBase


# Used to preserve row order in outputs.
ROW_ORDER_COLUMN = '_sqlite_id'


class TabularFileLoaderSQLite(TabularFileLoaderBase):
    """Class for loading tabular files into a SQLite database.

    Extends the base class with additional SQLite-specific attributes - see __init__.
    """
    def __init__(self, file:str, connection:sqlite3.Connection, table_name:str, header=True, names=[]):
        super().__init__(file, header, names)

        self.connection = connection
        "sqlite3 Connection object."

        self.table_name = table_name
        "Table name used to load tabular file contents into the database."

        self.columns = None
        "Column names to be used in the table."

    def load(self):
        # Insert ROW_ORDER_COLUMN as first column with dummy value of 1 so the schema matches what's given by self._iter_indexed_rows().
        self.df_head.insert(0, ROW_ORDER_COLUMN, value=1)
        self.columns = self.df_head.columns

        # Create table with schema defined by self.df_head which includes the additional ROW_ORDER_COLUMN.
        with self.connection:
            create_table_statement = pd.io.sql.get_schema(self.df_head, self.table_name, con=self.connection)
            self.connection.execute(create_table_statement)

        insert_statement = f"""
            INSERT INTO {self.table_name}
            VALUES ({','.join(['?' for _ in self.columns])})
        """
        rows = self._iter_indexed_rows()
        try:
            with self.connection:
                self.connection.executemany(insert_statement, rows)
        except sqlite3.ProgrammingError as e:
            raise ValueError(f'Failed to load {self.file}.') from e


def cleanup(database:str):
    """Removes the database file if present."""
    try:
        os.remove(database)
    except FileNotFoundError:
        pass


def sanitize_identifier(identifier:str):
    """Sanitize a SQLite identifier.
    
    1. Escape existing double quotes
    2. Wrap inside double quotes
    """
    identifier = identifier.replace('"', '""')
    return f'"{identifier}"'


def chunked_query_to_csv(con:sqlite3.Connection, query:str, path:str, chunksize:int, header=True, columns_to_exclude:List[str]=[], **to_csv_kwargs):
    """Write query results to a file in chunks so that query results do not need to be loaded into memory all at once.

    Parameters
    ----------
    con
        sqlite3 Connection object.
    query
        SQLite query string.
    path
        Path to output file.
    chunksize
        Size of chunks to load into memory.
    header
        Write out the column names.
    columns_to_exclude
        Column names to exclude from output.
        There is no way to exclude columns in a SQLite query, so it is implemented as an option here.
    to_csv_kwargs
        Additional keyword arguments passed to pandas.to_csv().
    """
    df_chunks = pd.read_sql_query(query, con, chunksize=chunksize)

    # Overwrite the file if it exists.
    mode = 'w'

    for df_chunk in df_chunks:
        if columns_to_exclude:
            for col in columns_to_exclude:
                df_chunk.drop(col, axis=1, inplace=True)
        df_chunk.to_csv(path, mode=mode, header=header, **to_csv_kwargs)

        # For subsequent chunks, use append mode and don't include the header.
        mode = 'a'
        header = False
