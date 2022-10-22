"""
A suite of commands to help with data curation.
"""
import argparse
import sys

from augur.argparse_ import add_command_subparsers
from augur.errors import AugurError
from augur.io.json import dump_ndjson, load_ndjson
from augur.io.metadata import read_table_to_dict
from augur.types import DataErrorMethod
from . import passthru


SUBCOMMAND_ATTRIBUTE = '_curate_subcommand'
SUBCOMMANDS = [
    passthru,
]


def create_shared_parser():
    """
    Creates an argparse.ArgumentParser that is intended to be used as a parent
    parser¹ for all `augur curate` subcommands. This should include all options
    that are intended to be shared across the subcommands.

    Note that any options strings used here cannot be used in individual subcommand
    subparsers unless the subparser specifically sets `conflict_handler='resolve'`²,
    then the subparser option will override the option defined here.

    Based on https://stackoverflow.com/questions/23296695/permit-argparse-global-flags-after-subcommand/23296874#23296874

    ¹ https://docs.python.org/3/library/argparse.html#parents
    ² https://docs.python.org/3/library/argparse.html#conflict-handler
    """
    shared_parser = argparse.ArgumentParser(add_help=False)

    shared_inputs = shared_parser.add_argument_group(
        title="INPUTS",
        description="""
            Input options shared by all `augur curate` commands.
            If no input options are provided, commands will try to read NDJSON records from stdin.
        """)
    shared_inputs.add_argument("--metadata",
        help="Input metadata file, as CSV or TSV.")
    shared_inputs.add_argument("--id-column",
        help="Name of the metadata column that contains the record identifier for reporting duplicate records. "
             "Uses the first column of the metadata file if not provided.")
    shared_inputs.add_argument("--duplicate-reporting",
        choices=[ method.value for method in DataErrorMethod ],
        default=DataErrorMethod.ERROR_FIRST.value,
        help="How should duplicate records be reported.")

    return shared_parser


def register_parser(parent_subparsers):
    shared_parser = create_shared_parser()
    parser = parent_subparsers.add_parser("curate", help=__doc__)

    # Add print_help so we can run it when no subcommands are called
    parser.set_defaults(print_help = parser.print_help)

    # Add subparsers for subcommands
    subparsers = parser.add_subparsers(dest="subcommand", required=False)
    # Add the shared_parser to make it available for subcommands
    # to include in their own parser
    subparsers.shared_parser = shared_parser
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
    if args.metadata:
        records = read_table_to_dict(args.metadata, DataErrorMethod(args.duplicate_reporting), args.id_column)
    elif not sys.stdin.isatty():
        records = load_ndjson(sys.stdin)
    else:
        raise AugurError("No valid inputs were provided.")

    # Run subcommand to get modified records
    modified_records = getattr(args, SUBCOMMAND_ATTRIBUTE).run(args, records)

    # Output modified records
    dump_ndjson(modified_records)
