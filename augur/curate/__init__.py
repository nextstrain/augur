"""
A suite of commands to help with data curation.
"""
from . import format_dates, normalize_strings, passthru, titlecase, apply_geolocation_rules, apply_record_annotations, abbreviate_authors, parse_genbank_location, transform_strain_name, rename


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

    parser.subcommands = SUBCOMMANDS

    return parser
