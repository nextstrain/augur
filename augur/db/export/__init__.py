"""
Export data from a database.
"""
from augur.errors import AugurError
from augur.utils import first_line
from .sqlite3 import export as export_sqlite3
from ..defaults import DEFAULT_IMPORTED_METADATA_TABLE, DEFAULT_IMPORTED_SEQUENCES_TABLE


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("export", help=first_line(__doc__))
    parser.add_argument('--db', help="SQLite3 database file to export from.")
    parser.add_argument('--output-metadata', help="Metadata file to output.")
    parser.add_argument('--output-sequences', help="Sequences file (FASTA) to output.")
    parser.add_argument('--metadata-table', default=DEFAULT_IMPORTED_METADATA_TABLE, help="Table name in database that contains metadata.")
    parser.add_argument('--sequences-table', default=DEFAULT_IMPORTED_SEQUENCES_TABLE, help="Table name in database that contains sequences.")
    return parser


def run(args):
    validate_args(args)

    export_sqlite3(args.db,
        args.metadata_table, args.sequences_table,
        args.output_metadata, args.output_sequences)


def validate_args(args):
    if not args.db:
        raise AugurError("Must provide a database filepath.")

    if not any((args.output_metadata, args.output_sequences)):
        raise AugurError("Must choose to output metadata and/or sequences.")
