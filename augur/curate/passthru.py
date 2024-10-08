"""
Pass through records without doing any data transformations.
Useful for testing, troubleshooting, or just converting file formats.
"""
from .argparse_shared_parser import shared_parser


def register_parser(parent_subparsers):
    return parent_subparsers.add_parser("passthru",
        parents=[shared_parser],
        help=__doc__)


def run(args, records):
    yield from records

