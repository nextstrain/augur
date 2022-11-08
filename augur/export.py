"""
Export JSON files suitable for visualization with auspice.
"""
from .argparse_ import add_command_subparsers
from . import export_v1, export_v2

SUBCOMMANDS = [
    export_v2,
    export_v1,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("export", help=__doc__)
    # Add subparsers for subcommands
    metavar_msg ="Augur export now needs you to define the JSON version " + \
                 "you want, e.g. `augur export v2`."
    subparsers = parser.add_subparsers(title="JSON SCHEMA",
                                       metavar=metavar_msg)
    add_command_subparsers(subparsers, SUBCOMMANDS)
    return parser
