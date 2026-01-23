"""
Concatenate multiple measurements JSONs into a single JSON file
"""
import sys

from augur.argparse_ import ExtendOverwriteDefault, add_minify_arguments
from augur.io.json import write_json
from augur.utils import first_line
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

    add_minify_arguments(parser)

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

    if args.minify_json:
        minify = True
    elif args.no_minify_json:
        minify = False
    else:
        minify = None
    write_json(output, args.output_json, minify=minify)
    try:
        read_measurements_json(measurements_json=args.output_json)
    except ValidateError:
        print(
            "ERROR: Validation of output JSON failed. See detailed errors above.",
            file=sys.stderr,
        )
        sys.exit(1)
