"""
Pass through records without doing any data transformations.
Useful for testing, troubleshooting, or just converting file formats.
"""


def register_parser(parent_subparsers, **kwargs):
    return parent_subparsers.add_parser("passthru",
        parents=[parent_subparsers.shared_parser],
        **kwargs)


def run(args, records):
    yield from records

