"""
Create JSON files suitable for visualization within the measurements panel of Auspice.
"""
from augur.argparse_ import add_command_subparsers
from . import export, concat

SUBCOMMANDS = [
    export,
    concat,
]


def register_parser(parent_subparsers, **kwargs):
    parser = parent_subparsers.add_parser("measurements", **kwargs)
    # Add subparsers for subcommands
    subparsers = parser.add_subparsers(dest='subcommand')
    add_command_subparsers(subparsers, SUBCOMMANDS)
    return parser
