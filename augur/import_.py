"""
Import analyses into augur pipeline from other systems
"""
from .utils import first_line
from .import_beast import run_beast, register_arguments_beast

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("import", help=first_line(__doc__))
    metavar_msg = "Import analyses into augur pipeline from other systems"
    subparsers = parser.add_subparsers(title="TYPE",
                                       metavar=metavar_msg)
    subparsers.required = True
    register_arguments_beast(subparsers)
    return parser

def run(args):
    if "beast" in args:
        return run_beast(args)
