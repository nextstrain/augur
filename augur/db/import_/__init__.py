"""
Import data into a database.
"""
from augur.errors import AugurError
from augur.utils import first_line
from .sqlite3 import import_ as import_sqlite3
from ..defaults import DEFAULT_IMPORTED_METADATA_TABLE, DEFAULT_IMPORTED_SEQUENCES_TABLE


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("import", help=first_line(__doc__))
    parser.add_argument('--metadata', '-m', help="Metadata file to import into the database.")
    parser.add_argument('--sequences', '-s', help="Sequences file (FASTA) to import into the database.")
    parser.add_argument('--output-db', help="Destination SQLite3 database file.")
    parser.add_argument('--metadata-table', default=DEFAULT_IMPORTED_METADATA_TABLE, help="Table name to store metadata in.")
    parser.add_argument('--sequences-table', default=DEFAULT_IMPORTED_SEQUENCES_TABLE, help="Table name to store sequences in.")
    return parser


def run(args):
    validate_args(args)

    import_sqlite3(args.metadata, args.metadata_table,
        args.sequences, args.sequences_table,
        args.output_db)


def validate_args(args):
    if not args.output_db:
        raise AugurError("Must provide a database filepath.")

    if not any((args.metadata, args.sequences)):
        raise AugurError("Must provide metadata and/or sequences.")
