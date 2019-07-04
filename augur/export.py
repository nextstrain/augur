"""
Export JSON files suitable for visualization with auspice.
"""
from .export_v1 import run_v1, register_arguments_v1
from .export_v2 import run_v2, register_arguments_v2


def register_arguments(parser):
    metavar_msg ="Augur export now needs you to define the JSON version " + \
                 "you want, e.g. `augur export v2`."
    subparsers = parser.add_subparsers(title="JSON SCHEMA",
                                       metavar=metavar_msg)
    subparsers.required = True
    register_arguments_v2(subparsers)
    register_arguments_v1(subparsers)


def run(args):
    if "v1" in args:
        return run_v1(args)
    else:
        return run_v2(args)
