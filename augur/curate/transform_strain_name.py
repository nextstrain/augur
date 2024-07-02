#!/usr/bin/env python3
"""
Verifies strain name pattern in the 'strain' field of the NDJSON record from
stdin. Adds a 'strain' field to the record if it does not already exist.

Outputs the modified records to stdout.
"""
import argparse
import json
import re
from sys import stderr, stdin, stdout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--strain-regex", default="^.+$",
        help="Regex pattern for strain names. " +
             "Strain names that do not match the pattern will be dropped.")
    parser.add_argument("--backup-fields", nargs="*",
        help="List of backup fields to use as strain name if the value in 'strain' " +
             "does not match the strain regex pattern. " +
             "If multiple fields are provided, will use the first field that has a non-empty string.")

    args = parser.parse_args()

    strain_name_pattern = re.compile(args.strain_regex)

    for index, record in enumerate(stdin):
        record = json.loads(record)

        # Verify strain name matches the strain regex pattern
        if strain_name_pattern.match(record.get('strain', '')) is None:
            # Default to empty string if not matching pattern
            record['strain'] = ''
            # Use non-empty value of backup fields if provided
            if args.backup_fields:
                for field in args.backup_fields:
                    if record.get(field):
                        record['strain'] = str(record[field])
                        break

        if record['strain'] == '':
            print(f"WARNING: Record number {index} has an empty string as the strain name.", file=stderr)


        json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
        print()
