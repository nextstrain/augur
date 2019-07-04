"""
Validate a set of JSON files intended for visualization in auspice.
"""
from .validate_v1 import run_v1, register_arguments_v1
from .validate_v2 import run_v2, register_arguments_v2

def register_arguments(parser):
    metavar_msg ="Augur validate now needs you to define the JSON version " + \
                 "you want to validate, e.g. `augur validate v2`."
    subparsers = parser.add_subparsers(title="JSON SCHEMA",
                                       metavar=metavar_msg)
    subparsers.required = True
    register_arguments_v1(subparsers)
    register_arguments_v2(subparsers)

def run(args):
    '''
    Validate auspice-compatable JSONs against a schema
    '''
    if "v1" in args:
        run_v1(args)
    else:
        run_v2(args)
