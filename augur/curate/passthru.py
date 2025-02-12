"""
Pass through records without doing any data transformations.
Useful for testing, troubleshooting, or just converting file formats.
"""
from augur.utils import first_line


def register_parser(parent_subparsers):
    return parent_subparsers.add_parser("passthru",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))


def run(args, records):
    yield from records
