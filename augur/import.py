"""
Import analyses into augur pipeline from other systems
"""
from .import_beast import run_beast, register_arguments_beast

def register_arguments(parser):
    metavar_msg = "Import analyses into augur pipeline from other systems"
    subparsers = parser.add_subparsers(title="TYPE",
                                       metavar=metavar_msg)
    subparsers.required = True
    register_arguments_beast(subparsers)

def run(args):
    if "beast" in args:
        return run_beast(args)
