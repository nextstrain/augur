"""
A suite of commands to help with data curation.
"""
import sys

from augur.argparse_ import add_command_subparsers
from augur.errors import AugurError
from augur.io.json import dump_ndjson, load_ndjson
from . import passthru


SUBCOMMAND_ATTRIBUTE = '_curate_subcommand'
SUBCOMMANDS = [
    passthru,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("curate",
        help=__doc__)

    # Add print_help so we can run it when no subcommands are called
    parser.set_defaults(print_help = parser.print_help)

    # Add subparsers for subcommands
    subparsers = parser.add_subparsers(dest="subcommand", required=False)
    # Using a subcommand attribute so subcommands are not directly
    # run by top level Augur. Process I/O in `curate`` so individual
    # subcommands do not have to worry about it.
    add_command_subparsers(subparsers, SUBCOMMANDS, SUBCOMMAND_ATTRIBUTE)

    return parser


def run(args):
    # Print help if no subcommands are used
    if not getattr(args, SUBCOMMAND_ATTRIBUTE, None):
        return args.print_help()

    # Read inputs
    if not sys.stdin.isatty():
        records = load_ndjson(sys.stdin)
    else:
        raise AugurError("No valid inputs were provided.")

    # Run subcommand to get modified records
    modified_records = getattr(args, SUBCOMMAND_ATTRIBUTE).run(args, records)

    # Output modified records
    dump_ndjson(modified_records)
