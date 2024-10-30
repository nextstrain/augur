"""
Import analyses into augur pipeline from other systems
"""
from augur.argparse_ import add_command_subparsers

SUBCOMMANDS = [
    "import_.beast",
]

def register_parser(parser):
    metavar_msg = "Import analyses into augur pipeline from other systems"
    subparsers = parser.add_subparsers(title="TYPE",
                                       metavar=metavar_msg)
    add_command_subparsers(subparsers, SUBCOMMANDS)
    return parser
