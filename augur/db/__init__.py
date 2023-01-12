"""
Commands to interface with a database.
"""
from augur.argparse_ import add_command_subparsers
from augur.utils import first_line
from . import import_, export

SUBCOMMANDS = [
    import_,
    export,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("db", help=first_line(__doc__))
    # Add subparsers for subcommands
    subparsers = parser.add_subparsers(dest='subcommand')
    add_command_subparsers(subparsers, SUBCOMMANDS)
    return parser
