"""
Applies record annotations to overwrite field values.
This does not do any additional transformations on top of the annotations.
"""
import csv
from collections import defaultdict
from augur.errors import AugurError
from augur.io.print import print_err
from augur.utils import first_line


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("apply-record-annotations",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))

    parser.add_argument("--annotations", metavar="TSV", required=True,
        help="Manually curated annotations TSV file. " +
             "The TSV should not have a header and should have exactly three columns: " +
             "id to match existing metadata, field name, and field value. " +
             "If there are multiple annotations for the same id and field, then the last value is used. " +
             "Lines starting with '#' are treated as comments. " +
             "Any '#' after the field value are treated as comments.")
    parser.add_argument("--id-field", default="accession",
        help="The ID field in the metadata to use to merge with the annotations.")

    return parser


def run(args, records):
    annotations = defaultdict(dict)
    with open(args.annotations, 'r', newline='') as annotations_fh:
        csv_reader = csv.reader(annotations_fh, delimiter='\t')
        for row in csv_reader:
            if not row or row[0].lstrip()[0] == '#':
                    continue
            elif len(row) != 3:
                print_err("WARNING: Could not decode annotation line " + "\t".join(row))
                continue
            id, field, value = row
            annotations[id][field] = value.partition('#')[0].rstrip()

    for record in records:
        record_id = record.get(args.id_field)
        if record_id is None:
            raise AugurError(f"ID field {args.id_field!r} does not exist in record")

        record_annotations = annotations.get(record_id, {})
        for field in list(record_annotations.keys()):
            if field not in record:
                print_err(f"WARNING: Skipping annotation for field {field!r} that does not exist in record")
                del record_annotations[field]

        record.update(record_annotations)

        yield record
