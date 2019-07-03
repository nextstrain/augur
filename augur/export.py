"""
Export JSON files suitable for visualization with auspice.
"""
from .export_v1 import run_v1, register_arguments_v1
from .export_v2 import run_v2, register_arguments_v2


def register_arguments(parser):
    subparsers = parser.add_subparsers()
    register_arguments_v2(subparsers)
    register_arguments_v1(subparsers)


def run(args):
    if "v1" in args:
        return run_v1(args)
    else:
        return run_v2(args)
