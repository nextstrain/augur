"""
A suite of commands to help with data curation.
"""
from augur.argparse_ import add_command_subparsers
from . import format_dates, normalize_strings, passthru, titlecase, apply_geolocation_rules, apply_record_annotations, abbreviate_authors, parse_genbank_location, transform_strain_name, rename


SUBCOMMAND_ATTRIBUTE = '_curate_subcommand'
SUBCOMMANDS = [
    passthru,
    normalize_strings,
    format_dates,
    titlecase,
    apply_geolocation_rules,
    apply_record_annotations,
    abbreviate_authors,
    parse_genbank_location,
    transform_strain_name,
    rename,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("curate", help=__doc__)

    # Add print_help so we can run it when no subcommands are called
    parser.set_defaults(print_help = parser.print_help)

    # Add subparsers for subcommands
    subparsers = parser.add_subparsers(dest="subcommand", required=False)
    # Using a subcommand attribute so subcommands are not directly
    # run by top level Augur. Process I/O in `curate`` so individual
    # subcommands do not have to worry about it.
    add_command_subparsers(subparsers, SUBCOMMANDS, SUBCOMMAND_ATTRIBUTE)

    return parser
