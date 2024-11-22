"""
Create JSON files suitable for visualization within the measurements panel of Auspice.
"""
from augur.utils import first_line
from . import export, concat

SUBCOMMANDS = [
    export,
    concat,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("measurements", help=first_line(__doc__))

    parser.subcommands = SUBCOMMANDS

    return parser
