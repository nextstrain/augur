"""
Renames fields of the NDJSON record from stdin and outputs modified records
to stdout.
"""
import argparse
import json
from sys import stderr, stdin, stdout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--field-map", nargs="+",
        help="Fields names in the NDJSON record mapped to new field names, " +
             "formatted as '{old_field_name}={new_field_name}'. " +
             "If the old field does not exist in record, the new field will be added with an empty string value. " +
             "If the new field already exists in record, then the renaming of the old field will be skipped. " +
             "Skips the field if the old field name is the same as the new field name (case-sensitive).")
    parser.add_argument("--force", action="store_true",
        help="Force renaming of old field even if the new field already exists. " +
             "Please keep in mind this will overwrite the value of the new field.")

    args = parser.parse_args()

    field_map = {}
    for field in args.field_map:
        old_name, new_name = field.split('=')

        if old_name == new_name:
            continue

        field_map[old_name] = new_name

    for record in stdin:
        record = json.loads(record)

        for old_field, new_field in field_map.items():

            if record.get(new_field) and not args.force:
                print(
                    f"WARNING: skipping rename of {old_field} because record",
                    f"already has a field named {new_field}.",
                    file=stderr
                )
                continue

            record[new_field] = record.pop(old_field, '')

        json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
        print()
