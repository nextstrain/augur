"""
Import analyses into augur pipeline from other systems
"""
from augur.utils import first_line
from . import beast

SUBCOMMANDS = [
    beast,
]

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("import", help=first_line(__doc__))

    parser.subcommands = SUBCOMMANDS

    return parser
