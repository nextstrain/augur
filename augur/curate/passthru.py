"""
Pass through records without doing any data transformations.
Useful for testing, troubleshooting, or just converting file formats.
"""
from ._shared import shared_parser, validate


def register_parser(parent_subparsers):
    return parent_subparsers.add_parser("passthru",
        parents=[shared_parser],
        help=__doc__)


@validate
def run(args, records):
    yield from records

