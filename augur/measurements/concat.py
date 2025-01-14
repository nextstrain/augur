"""
Concatenate multiple measurements JSONs into a single JSON file
"""
import sys

from augur.argparse_ import ExtendOverwriteDefault
from augur.utils import first_line, write_json
from augur.validate import (
    measurements as read_measurements_json,
    ValidateError
)


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("concat", help=first_line(__doc__))

    required = parser.add_argument_group(
        title="REQUIRED"
    )
    required.add_argument("--jsons", required=True, type=str, nargs="+", action=ExtendOverwriteDefault, metavar="JSONs",
        help="Measurement JSON files to concatenate.")
    required.add_argument("--output-json", required=True, metavar="JSON", type=str,
        help="Output JSON file")

    optional = parser.add_argument_group(
        title="OPTIONAL SETTINGS"
    )
    optional.add_argument("--default-collection", type=str,
        help="The key of the default collection to display. " +
             "If not provided, the first collection of the first JSON file will be displayed")
    optional.add_argument("--minify-json", action="store_true",
        help="Concatenate JSONs without indentation or line returns.")

    return parser


def run(args):
    output = {
        'collections': []
    }
    if args.default_collection is not None:
        output['default_collection'] = args.default_collection

    for json in args.jsons:
        measurements = read_measurements_json(json)
        output['collections'].extend(measurements['collections'])

    indent = {"indent": None} if args.minify_json else {}
    write_json(output, args.output_json, include_version=False, **indent)
    try:
        read_measurements_json(measurements_json=args.output_json)
    except ValidateError:
        print(
            "ERROR: Validation of output JSON failed. See detailed errors above.",
            file=sys.stderr,
        )
        sys.exit(1)
