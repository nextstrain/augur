#!/usr/bin/env python3
"""
Merges user curated annotations with the NDJSON records from stdin, with the user
curations overwriting the existing fields. The modified records are output
to stdout. This does not do any additional transformations on top of the user
curations.
"""
import argparse
import csv
import json
from collections import defaultdict
from sys import exit, stdin, stderr, stdout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--annotations", metavar="TSV", required=True,
        help="Manually curated annotations TSV file. " +
             "The TSV should not have a header and should have exactly three columns: " +
             "id to match existing metadata, field name, and field value. " +
             "If there are multiple annotations for the same id and field, then the last value is used. " +
             "Lines starting with '#' are treated as comments. " +
             "Any '#' after the field value are treated as comments.")
    parser.add_argument("--id-field", default="accession",
        help="The ID field in the metadata to use to merge with the annotations.")

    args = parser.parse_args()

    annotations = defaultdict(dict)
    with open(args.annotations, 'r') as annotations_fh:
        csv_reader = csv.reader(annotations_fh, delimiter='\t')
        for row in csv_reader:
            if not row or row[0].lstrip()[0] == '#':
                    continue
            elif len(row) != 3:
                print("WARNING: Could not decode annotation line " + "\t".join(row), file=stderr)
                continue
            id, field, value = row
            annotations[id][field] = value.partition('#')[0].rstrip()

    for record in stdin:
        record = json.loads(record)

        record_id = record.get(args.id_field)
        if record_id is None:
            print(f"ERROR: ID field {args.id_field!r} does not exist in record", file=stderr)
            exit(1)

        record.update(annotations.get(record_id, {}))

        json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
        print()
