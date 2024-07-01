"""
Renames fields / columns of the input data
"""

from typing import Iterable
from augur.io.print import print_err
import argparse

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("rename",
    parents = [parent_subparsers.shared_parser],
    help = __doc__)

    p = parser.add_argument_group(title="RENAME SPECIFIC OPTIONS")

    p.add_argument("--field-map", nargs="+", default=[],
        help="Fields names in the NDJSON record mapped to new field names, " +
             "formatted as '{old_field_name}={new_field_name}'. " +
             "If the old field does not exist in record, the new field will be added with an empty string value. " +
             "If the new field already exists in record, then the renaming of the old field will be skipped. " +
             "Skips the field if the old field name is the same as the new field name (case-sensitive).")
    p.add_argument("--force", action="store_true",
        help="Force renaming of old field even if the new field already exists. " +
             "Please keep in mind this will overwrite the value of the new field.")

    return parser

def run(args: argparse.Namespace, records: Iterable[dict]) -> Iterable[dict]:

    field_map = {}
    for field in args.field_map:
        old_name, new_name = field.split('=')

        if old_name == new_name:
            continue

        field_map[old_name] = new_name

    for record in records:
        record = record.copy()

        for old_field, new_field in field_map.items():

            if record.get(new_field) and not args.force:
                print_err(
                    f"WARNING: skipping rename of {old_field} because record",
                    f"already has a field named {new_field}."
                )
                continue

            record[new_field] = record.pop(old_field, '')

        yield(record)
