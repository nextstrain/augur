import os
from augur.io.defaults import POTENTIAL_STRAIN_ID_COLUMNS

from augur.io.metadata import Metadata
from augur.io.sequences import Sequences
from augur.io.sqlite3 import Sqlite3Database, Sqlite3DatabaseError
from ..defaults import DEFAULT_SEQUENCES_ID_COLUMN, DEFAULT_SEQUENCES_SEQ_COLUMN


def import_(metadata_file: str, metadata_table_name: str,
        sequences_file: str, sequences_table_name: str,
        db_file: str):
    """Import data into a SQLite3 database file."""
    if metadata_file:
        db_file_existed = os.path.exists(db_file)
        with Sqlite3Database(db_file) as database:
            try:
                metadata = Metadata(metadata_file)
                import_metadata(metadata, metadata_table_name, database)
            except Sqlite3DatabaseError as e:
                # Delete the database file if it was created for this import.
                if not db_file_existed:
                    database.delete()
                raise e
    if sequences_file:
        db_file_existed = os.path.exists(db_file)
        with Sqlite3Database(db_file) as database:
            try:
                # Get ID column name from metadata if available.
                id_column = existing_id_column(database, metadata_table_name) or DEFAULT_SEQUENCES_ID_COLUMN

                sequences = Sequences(sequences_file)
                import_sequences(sequences, id_column, sequences_table_name, database)
            except Sqlite3DatabaseError as e:
                # Delete the database file if it was created for this import.
                if not db_file_existed:
                    database.delete()
                raise e


def import_metadata(metadata: Metadata, table_name: str, database: Sqlite3Database):
    """Import metadata into a SQLite3 database."""
    database.create_table(table_name, metadata.names)
    database.insert(table_name, metadata.names, metadata.rows())


def import_sequences(sequences: Sequences, id_column: str, table_name: str, database: Sqlite3Database):
    """Import sequences into a SQLite3 database."""
    columns = [id_column, DEFAULT_SEQUENCES_SEQ_COLUMN]
    database.create_table(table_name, columns)
    records = (
        {id_column: record.id, DEFAULT_SEQUENCES_SEQ_COLUMN: str(record.seq)}
        for record in sequences
    )
    database.insert(table_name, columns, records)


def existing_id_column(database: Sqlite3Database, table_name: str):
    existing_columns = set(database.columns(table_name))
    for col in POTENTIAL_STRAIN_ID_COLUMNS:
        if col in existing_columns:
            return col
