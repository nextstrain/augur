"""
Verifies strain name pattern in the 'strain' field.
Adds a 'strain' field to the record if it does not already exist.
"""

import argparse
import re
from typing import Generator, List
from augur.argparse_ import ExtendOverwriteDefault
from augur.io.print import print_err
from augur.utils import first_line


def transform_name(
    record: dict,
    index: int,
    strain_name_pattern: re.Pattern,
    backup_fields: List[str],
) -> dict:
    # Verify strain name matches the strain regex pattern
    if strain_name_pattern.match(record.get("strain", "")) is None:
        # Default to empty string if not matching pattern
        record["strain"] = ""

        # Use non-empty value of backup fields if provided
        if backup_fields:
            for field in backup_fields:
                if record.get(field):
                    record["strain"] = str(record[field])
                    break

    if record["strain"] == "":
        print_err(
            f"WARNING: Record number {index} has an empty string as the strain name.",
        )

    return record


def register_parser(
    parent_subparsers: argparse._SubParsersAction,
) -> argparse._SubParsersAction:
    parser = parent_subparsers.add_parser(
        "transform-strain-name",
        parents=[parent_subparsers.shared_parser],  # type: ignore[attr-defined]
        help=first_line(__doc__),
    )

    parser.add_argument(
        "--strain-regex",
        default="^.+$",
        help="Regex pattern for strain names. "
        + "Strain names that do not match the pattern will be dropped.",
    )
    parser.add_argument(
        "--backup-fields",
        nargs="*",
        action=ExtendOverwriteDefault,
        default=[],
        help="List of backup fields to use as strain name if the value in 'strain' "
        + "does not match the strain regex pattern. "
        + "If multiple fields are provided, will use the first field that has a non-empty string.",
    )

    return parser


def run(args: argparse.Namespace, records: List[dict]) -> Generator[dict, None, None]:
    strain_name_pattern = re.compile(args.strain_regex)

    for index, record in enumerate(records):
        transform_name(
            record,
            index,
            strain_name_pattern,
            args.backup_fields,
        )

        yield record
