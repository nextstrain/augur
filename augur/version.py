"""
Print the version of augur.
"""
from .__version__ import __version__

def register_parser(parent_subparsers, **kwargs):
    return parent_subparsers.add_parser("version", **kwargs)

def run(args):
    print("augur", __version__)
    return 0
