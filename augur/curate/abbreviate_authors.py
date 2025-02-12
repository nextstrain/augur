"""
Abbreviates a full list of authors to be '<first author> et al.'
Expects NDJSON records from stdin and outputs modified records to stdout.

Note: This is a "best effort" approach and can potentially mangle the author name.
"""

import argparse
import re
from typing import Generator, List
from augur.io.print import print_err
from augur.utils import first_line


def parse_authors(
    record: dict,
    authors_field: str,
    default_value: str,
    index: int,
    abbr_authors_field: str = None,
) -> dict:
    # Strip and normalize whitespace
    new_authors = re.sub(r"\s+", " ", record[authors_field])

    if new_authors == "":
        new_authors = default_value
    else:
        # Split authors list on comma/semicolon
        # OR "and"/"&" with at least one space before and after
        new_authors = re.split(r"(?:\s*[,，;；]\s*|\s+(?:and|&)\s+)", new_authors)[0]

        # if it does not already end with " et al.", add it
        if not new_authors.strip(". ").endswith(" et al"):
            new_authors += " et al."

    if abbr_authors_field:
        if record.get(abbr_authors_field):
            print_err(
                f"WARNING: the {abbr_authors_field!r} field already exists",
                f"in record {index} and will be overwritten!",
            )

        record[abbr_authors_field] = new_authors
    else:
        record[authors_field] = new_authors

    return record


def register_parser(
    parent_subparsers: argparse._SubParsersAction,
) -> argparse._SubParsersAction:
    parser = parent_subparsers.add_parser(
        "abbreviate-authors",
        parents=[parent_subparsers.shared_parser],  # type: ignore[attr-defined]
        help=first_line(__doc__),
    )

    parser.add_argument(
        "--authors-field",
        default="authors",
        help="The field containing list of authors.",
    )
    parser.add_argument(
        "--default-value",
        default="?",
        help="Default value to use if authors list is empty.",
    )
    parser.add_argument(
        "--abbr-authors-field",
        help="The field for the generated abbreviated authors. "
        + "If not provided, the original authors field will be modified.",
    )

    return parser


def run(args: argparse.Namespace, records: List[dict]) -> Generator[dict, None, None]:
    for index, record in enumerate(records):
        parse_authors(
            record,
            args.authors_field,
            args.default_value,
            index,
            args.abbr_authors_field,
        )

        yield record
