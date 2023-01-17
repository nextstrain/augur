"""
Import analyses into augur pipeline from other systems
"""
from augur.argparse_ import add_command_subparsers
from . import beast

SUBCOMMANDS = [
    beast,
]

def register_parser(parent_subparsers, **kwargs):
    parser = parent_subparsers.add_parser("import", **kwargs)
    metavar_msg = "Import analyses into augur pipeline from other systems"
    subparsers = parser.add_subparsers(title="TYPE",
                                       metavar=metavar_msg)
    add_command_subparsers(subparsers, SUBCOMMANDS)
    return parser
