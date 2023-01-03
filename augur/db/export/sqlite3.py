from augur.io.sqlite3 import Sqlite3Database
from ..defaults import DEFAULT_SEQUENCES_ID_COLUMN, DEFAULT_SEQUENCES_SEQ_COLUMN


def export(db_file: str,
        metadata_table: str, sequences_table: str,
        output_metadata: str = None, output_sequences: str = None):
    """Export data from a SQLite3 database file."""
    database = Sqlite3Database(db_file)

    if output_metadata:
        export_metadata(database, metadata_table, output_metadata)

    if output_sequences:
        export_sequences(database, sequences_table, output_sequences)


def export_metadata(database: Sqlite3Database, table_name: str, output_file: str):
    """Export metadata from a SQLite3 database."""

    database.table_to_csv(
        table_name=table_name,
        path=output_file,
        header=True,
        index=False,
        sep='\t',
    )


def export_sequences(database: Sqlite3Database, table_name: str, output_file: str):
    """Export sequences from a SQLite3 database."""

    # TODO: detect column names
    database.table_to_fasta(
        table_name=table_name,
        path=output_file,
        id_column=DEFAULT_SEQUENCES_ID_COLUMN,
        sequence_column=DEFAULT_SEQUENCES_SEQ_COLUMN,
    )
