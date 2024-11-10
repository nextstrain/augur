"""
Export a measurements JSON for a single collection
"""
import os
import pandas as pd
import sys

from augur.argparse_ import HideAsFalseAction
from augur.io.file import PANDAS_READ_CSV_OPTIONS
from augur.utils import first_line, write_json
from augur.validate import (
    measurements as read_measurements_json,
    measurements_collection_config as read_collection_config_json,
    ValidateError
)

# Default values for optional arguments that can also be provided via config file
# Setting as global dict instead of using argparse default so that the
# config file does not always get overwritten by the default values
DEFAULT_ARGS = {
    'title': 'Measurements',
    'x_axis_label': 'measurement values',
}


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("export", help=first_line(__doc__))

    required = parser.add_argument_group(
        title="REQUIRED"
    )
    required.add_argument("--collection", required=True, metavar="TSV",
        help="Collection of measurements and metadata in a TSV file. " +
             "Keep in mind duplicate columns will be renamed as 'X', 'X.1', 'X.2'...'X.N'")
    required.add_argument("--strain-column", default="strain",
        help="Name of the column containing strain names. " +
             "Provided column will be renamed to `strain` so please make sure no other columns are named `strain`. " +
             "Strain names in this column should match the strain names in the corresponding Auspice dataset JSON. " +
             "(default: %(default)s)")
    required.add_argument("--value-column", default="value",
        help="Name of the column containing the numeric values to be plotted for the given collection. " +
             "Provided column will be renamed to `value` so please make sure no other columns are named `value`. " +
             "(default: %(default)s)")
    required.add_argument("--output-json", required=True, metavar="JSON", type=str,
        help="Output JSON file. " +
             "The file name must follow the Auspice sidecar file naming convention to be recognized as a sidecar file. " +
             "See Nextstrain data format docs for more details.")

    config = parser.add_argument_group(
        title="COLLECTION CONFIGURATION",
        description="These options control the configuration of the collection for Auspice. " +
                    "You can provide a config JSON (which includes all available options) or " +
                    "command line arguments (which are more limited). " +
                    "Command line arguments will override the values set in the config JSON."
    )
    config.add_argument("--collection-config", metavar="JSON",
        help="Collection configuration file for advanced configurations. ")
    config.add_argument("--grouping-column", nargs="+", action="extend",
        help="Name of the column(s) that should be used as grouping(s) for measurements. " +
             "Note that if groupings are provided via command line args, the default group-by " +
             "field in the config JSON will be dropped.")
    config.add_argument("--key",
        help="A short key name of the collection for internal use within Auspice. " +
             "If not provided via config or command line option, the collection TSV filename will be used. ")
    config.add_argument("--title",
        help="The full title of the collection to display in the measurements panel title. " +
             f"If not provided via config or command line option, the panel's default title is {DEFAULT_ARGS['title']!r}.")
    config.add_argument("--x-axis-label",
        help="The short label to display for the x-axis that describles the value of the measurements. " +
             "If not provided via config or command line option, the panel's default " +
             f"x-axis label is {DEFAULT_ARGS['x_axis_label']!r}.")
    config.add_argument("--thresholds", type=float, nargs="+", action="extend",
        help="Measurements value threshold(s) to be displayed in the measurements panel.")
    config.add_argument("--filters", nargs="+", action="extend",
        help="The columns that are to be used a filters for measurements. " +
             "If not provided, all columns will be available as filters.")
    config.add_argument("--group-by", type=str,
        help="The default grouping column. If not provided, the first grouping will be used.")
    config.add_argument("--measurements-display", type=str, choices=["raw", "mean"],
        help="The default display of the measurements")

    config.add_argument("--show-overall-mean", "--hide-overall-mean",
        dest="show_overall_mean", action=HideAsFalseAction, nargs=0,
        help="Show or hide the overall mean per group by default")
    config.add_argument("--show-threshold", "--hide-threshold",
        dest="show_threshold", action=HideAsFalseAction, nargs=0,
        help="Show or hide the threshold(s) by default. This will be ignored if no threshold(s) are provided.")

    optional = parser.add_argument_group(
        title="OPTIONAL SETTINGS"
    )
    optional.add_argument("--include-columns", nargs="+", action="extend",
        help="The columns to include from the collection TSV in the measurements JSON. " +
             "Be sure to list columns that are used as groupings and/or filters. " +
             "If no columns are provided, then all columns will be included by default.")
    optional.add_argument("--minify-json", action="store_true",
        help="Export JSON without indentation or line returns.")

    return parser


def run(args):
    # Default value to None so all columns will be read
    columns_to_include = None
    if args.include_columns is not None:
        columns_to_include = set([args.strain_column, args.value_column] + args.include_columns)

    # Load input collection TSV file
    try:
        collection_df = pd.read_csv(args.collection, sep="\t", usecols=columns_to_include, **PANDAS_READ_CSV_OPTIONS)
    except FileNotFoundError:
        print(
            f"ERROR: collection TSV file {args.collection!r} does not exist",
            file=sys.stderr,
        )
        sys.exit(1)

    # Verify the strain and value columns are different
    if args.strain_column == args.value_column:
        print(
            "ERROR: The strain column and value column cannot be the same column.",
            file=sys.stderr
        )
        sys.exit(1)

    # Define mapping of required columns to user provided columns
    required_column_map = {
        args.strain_column: 'strain',
        args.value_column: 'value',
    }

    # Check all required columns are included in collection TSV
    checks_passed = True
    for provided_column, required_column in required_column_map.items():
        # Confirm the provided column exists
        if provided_column not in collection_df.columns:
            print(
                f"ERROR: Provided {required_column} column {provided_column!r} does not exist in collection TSV.",
                file=sys.stderr,
            )
            checks_passed = False
        # Confirm the provided column does not overwrite an existing column
        if (required_column in collection_df.columns and
            provided_column != required_column):
            print(
                f"ERROR: Cannot use provided {provided_column!r} column as the {required_column} column",
                f"because a {required_column!r} column already exists in collection TSV.",
                file=sys.stderr,
            )
            checks_passed = False

    if not checks_passed:
        sys.exit(1)

    # Rename user provided columns to expected columns
    collection_df = collection_df.rename(columns=required_column_map)

    # Make sure the value column is numeric
    try:
        collection_df['value'] = pd.to_numeric(collection_df['value'])
    except ValueError as e:
        print(f"ERROR: Found a non-numeric measurement value: {e!r}", file=sys.stderr)
        sys.exit(1)

    collection_config = {}
    if args.collection_config is not None:
        try:
            collection_config = read_collection_config_json(args.collection_config)
        except ValidateError:
            print(
                f"Validation of provided collection config JSON {args.collection_config!r} failed.",
                "Please check the formatting of this file.",
                file=sys.stderr
            )
            sys.exit(1)

    groupings = collection_config.pop('groupings', None)
    if args.grouping_column is not None:
        groupings = []
        for column in args.grouping_column:
            # If the user specified columns to include, verify the grouping column was included
            if args.include_columns and column not in args.include_columns:
                print(
                    f"ERROR: Provided grouping column {column!r} was not in the",
                    f"list of columns to include: {args.include_columns}.",
                    file=sys.stderr,
                )
                sys.exit(1)

            # Verify the grouping column is included in the collection
            if column not in collection_df.columns:
                print(
                    f"ERROR: Provided grouping column {column!r} does not exist in collection TSV.",
                    file=sys.stderr,
                )
                sys.exit(1)

            groupings.append({'key': column})

        if collection_config.get('display_defaults', {}).pop('group_by', None):
            print(
                "WARNING: The default group-by in the collection config has been removed",
                "because new groupings have been provided via the --grouping-column option.",
                file=sys.stderr
            )

    if not groupings:
        print("ERROR: Cannot create measurements JSON without valid groupings", file=sys.stderr)
        sys.exit(1)

    # Convert deprecated single threshold value to list if "thresholds" not provided
    single_threshold = collection_config.pop("threshold", None)
    if single_threshold is not None and "thresholds" not in collection_config:
        collection_config["thresholds"] = [single_threshold]

    # Combine collection config with command line args
    config_key_args = ['key', 'title', 'filters', 'x_axis_label', 'thresholds']
    display_default_args = ['group_by', 'measurements_display', 'show_overall_mean', 'show_threshold']

    for key_arg in config_key_args:
        key_arg_value = getattr(args, key_arg)
        if key_arg_value is not None:
            collection_config[key_arg] = key_arg_value

    for default_arg in display_default_args:
        default_arg_value = getattr(args, default_arg)
        if default_arg_value is not None:
            collection_config['display_defaults'] = collection_config.get('display_defaults', {})
            collection_config['display_defaults'][default_arg] = default_arg_value

    # Create collection output object with default values for required keys
    collection_output = {
        'key': collection_config.pop('key', os.path.basename(args.collection)),
        'title': collection_config.pop('title', DEFAULT_ARGS['title']),
        'groupings': groupings,
        'x_axis_label': collection_config.pop('x_axis_label', DEFAULT_ARGS['x_axis_label']),
        'measurements': collection_df.to_dict(orient='records'),
        **collection_config
    }

    # Create final output with single collection
    output = {
        'collections': [collection_output]
    }

    # Set indentation to None to create compact JSON if specified
    indent = {"indent": None} if args.minify_json else {}
    # Create output JSON
    write_json(output, args.output_json, include_version=False, **indent)
    # Verify the produced output is a valid measurements JSON
    try:
        read_measurements_json(measurements_json=args.output_json)
    except ValidateError:
        print(
            "ERROR: Validation of output JSON failed. See detailed errors above.",
            file=sys.stderr,
        )
        sys.exit(1)
