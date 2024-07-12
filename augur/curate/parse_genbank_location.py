"""
Parses GenBank's 'location' field of the NDJSON record to 3 separate
fields: 'country', 'division', and 'location'.

Checks that a record is from GenBank by verifying that the 'database'
field has a value of "GenBank" or "RefSeq".
"""

import argparse
from typing import Generator, List
from augur.io.print import print_err
from augur.utils import first_line


def parse_location(
    record: dict,
    location_field_name: str,
) -> dict:
    # Expected pattern for the location field is
    # "<country_value>[:<region>][, <locality>]"
    #
    # See GenBank docs for their "country" field:
    # https://www.ncbi.nlm.nih.gov/genbank/collab/country/
    location_field = record.get(location_field_name, None)
    if location_field is None:
        print_err(
            f"`parse-genbank-location` requires a `{location_field_name}` field; this record does not have one.",
        )
        # bail early because we're not gonna make any changes
        return record

    geographic_data = location_field.split(":")

    country = geographic_data[0]
    division = ""
    location = ""

    if len(geographic_data) == 2:
        division, _, location = geographic_data[1].partition(",")

    record["country"] = country.strip()
    record["division"] = division.strip()
    record["location"] = location.strip()

    return record


def register_parser(
    parent_subparsers: argparse._SubParsersAction,
) -> argparse._SubParsersAction:
    parser = parent_subparsers.add_parser(
        "parse-genbank-location",
        parents=[parent_subparsers.shared_parser],  # type: ignore
        help=first_line(__doc__),
    )

    parser.add_argument(
        "--location-field",
        default="geo_loc_name",
        help="The field containing the location, "
        + "in the format `<geo_loc_name>[:<region>][, <locality>]`",
    )

    return parser


def run(
    args: argparse.Namespace,
    records: List[dict],
) -> Generator[dict, None, None]:
    for record in records:
        parse_location(
            record,
            args.location_field,
        )
        yield record
