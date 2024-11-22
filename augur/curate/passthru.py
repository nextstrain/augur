"""
Pass through records without doing any data transformations.
Useful for testing, troubleshooting, or just converting file formats.
"""
from ._shared import shared_parser, validate


COMMAND_NAME = "passthru"


def register_parser(parent_subparsers):
    return parent_subparsers.add_parser(COMMAND_NAME,
        parents=[shared_parser],
        help=__doc__)


@validate(COMMAND_NAME)
def run(args, records):
    yield from records

