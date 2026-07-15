"""
Export JSON files suitable for visualization with auspice.
"""
from . import export_v1, export_v2

SUBCOMMANDS = [
    export_v2,
    export_v1,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("export", help=__doc__)
    parser.subcommands = SUBCOMMANDS
    return parser
