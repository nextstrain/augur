"""
Print the version of augur.
"""
from .__version__ import __version__

def register_parser(parent_subparsers):
    return parent_subparsers.add_parser("version", help=__doc__)

def run(args):
    print("augur", __version__)
    return 0
