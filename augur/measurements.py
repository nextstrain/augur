"""
Create JSON files suitable for visualization within the measurements panel of Auspice.
"""
import os
import pandas as pd
import sys

from .utils import write_json, HideAsFalseAction
from .validate import (
    measurements as validate_measurements_json,
    measurements_collection_config as validate_collection_config_json,
    ValidateError
)

# Default values for optional arguments that can also be provided via config file
# Setting as global dict instead of using argparse default so that the
# config file does not always get overwritten by the default values
DEFAULT_ARGS = {
    'title': 'Measurements',
    'x_axis_label': 'measurement values',
}


def load_collection(collection, strain_column, value_column):
    """
    Loads the provided collection TSV as a pandas DataFrame.
    Renames the provided strain and value columns if needed and ensures the
    value column has a numeric dtype.

    Prints any error messages to stderr.

    Parameters
    ----------
    collection: str
        Filepath to the collection TSV file
    strain_column: str
        The name of the strain column within the collection TSV
    value_column: str
        The name of the value column within the collection TSV

    Returns
    -------
    pandas.DataFrame or None
        The collection DataFrame or None if any errors were encountered during loading
    """
    try:
        collection_df = pd.read_csv(collection, sep="\t")
    except FileNotFoundError:
        print(
            f"ERROR: collection TSV file {collection!r} does not exist",
            file=sys.stderr,
        )
        return None

    # Verify the strain and value columns are different
    if strain_column == value_column:
        print(
            "ERROR: The strain column and value column cannot be the same column.",
            file=sys.stderr
        )
        return None

    # Define mapping of requried columns to user provided columns
    required_column_map = {
        strain_column: 'strain',
        value_column: 'value',
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
        return None

    # Rename user provided columns to expected columns
    collection_df = collection_df.rename(columns=required_column_map)

    # Make sure the value column is numeric
    try:
        collection_df['value'] = pd.to_numeric(collection_df['value'])
    except ValueError as e:
        print(f"ERROR: Found a non-numeric measurement value: {e!r}", file=sys.stderr)
        return None

    return collection_df


def get_collection_groupings(collection, grouping_columns):
    """
    Creates the groupings for the provided collection using the provided
    grouping columns after verifying the columns exist in the collection.

    Parameters
    ----------
    collection: pandas.DataFrame
        The collection used to validate groupings
    grouping_columns: list[str]
        List of grouping column names

    Returns
    -------
    list[dict] or None
        The groupings for the collection config or None any grouping columns are invalid
    """
    groupings = []
    for column in grouping_columns:
        if column not in collection.columns:
            print(
                f"ERROR: Provided grouping column {column!r} does not exist in collection TSV.",
                file=sys.stderr,
            )
            return None

        groupings.append({'key': column})

    return groupings


def override_config_with_args(config, args):
    """
    Overrides values in the config with values of provided command line args.

    Parameters
    ----------
    config: dict
        A collection config
    args: dict
        The __dict__ attribute of the parsed arguments from argparse
    """
    config_key_args = ['key', 'title', 'filters', 'x_axis_label', 'threshold']
    display_default_args = ['group_by', 'measurements_display', 'show_overall_mean', 'show_threshold']

    for key_arg in config_key_args:
        if args.get(key_arg) is not None:
            config[key_arg] = args[key_arg]

    for default_arg in display_default_args:
        if args.get(default_arg) is not None:
            config['display_defaults'] = config.get('display_defaults', {})
            config['display_defaults'][default_arg] = args[default_arg]


def validate_output_json(output_json):
    """
    Validate the output JSON against the measurements schema

    Parameters
    ----------
    output_json: str
        Filepath to output JSON

    """
    print("Validating produced measurements JSON")
    try:
        validate_measurements_json(measurements_json=output_json)
    except ValidateError:
        print(
            "ERROR: Validation of output JSON failed. See detailed errors above.",
            file=sys.stderr,
        )
        sys.exit(1)


def export_measurements(args):
    # Load input collection TSV file
    collection_df = load_collection(args['collection'], args['strain_column'], args['value_column'])

    if collection_df is None:
        print("ERROR: Loading of collection TSV was unsuccessful. See detailed errors above.", file=sys.stderr)
        sys.exit(1)

    collection_config = {}
    if args.get('collection_config'):
        try:
            collection_config = validate_collection_config_json(args['collection_config'])
        except ValidateError:
            print(
                f"Validation of provided collection config JSON {args['collection_config']!r} failed.",
                "Please check the formatting of this file.",
                file=sys.stderr
            )
            sys.exit(1)

    groupings = collection_config.pop('groupings', None)
    if args.get('grouping_column'):
        groupings = get_collection_groupings(collection_df, args['grouping_column'])
        if collection_config.get('display_defaults', {}).pop('group_by', None):
            print(
                "WARNING: The default group-by in the collection config has been removed",
                "because new groupings have been provided via the --grouping-column option.",
                file=sys.stderr
            )

    if not groupings:
        print("ERROR: Cannot create measurements JSON without valid groupings", file=sys.stderr)
        sys.exit(1)

    # Combine collection config with command line args
    override_config_with_args(collection_config, args)

    # Create collection output object with default values for required keys
    collection_output = {
        'key': collection_config.pop('key', os.path.basename(args['collection'])),
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
    indent = {"indent": None} if args['minify_json'] else {}
    # Create output JSON
    write_json(output, args['output_json'], include_version=False, **indent)
    # Verify the produced output is a valid measurements JSON
    validate_output_json(args['output_json'])


def concat_measurements(args):
    output = {
        'collections': []
    }
    if args.get("default_collection"):
        output['default_collection'] = args['default_collection']

    for json in args['jsons']:
        measurements = validate_measurements_json(json)
        output['collections'].extend(measurements['collections'])

    indent = {"indent": None} if args['minify_json'] else {}
    write_json(output, args['output_json'], include_version=False, **indent)
    validate_output_json(args['output_json'])


def register_arguments(parser):
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    export = subparsers.add_parser("export", help="Export a measurements JSON for a single collection")

    required = export.add_argument_group(
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

    config = export.add_argument_group(
        title="COLLECTION CONFIGURATION",
        description="These options control the configuration of the collection for Auspice. " +
                    "You can provide a config JSON (which includes all available options) or " +
                    "command line arguments (which are more limited). " +
                    "Command line arguments will override the values set in the config JSON."
    )
    config.add_argument("--collection-config", metavar="JSON",
        help="Collection configuration file for advanced configurations. ")
    config.add_argument("--grouping-column", nargs="+",
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
    config.add_argument("--threshold", type=float,
        help="A measurements value threshold to be displayed in the measurements panel.")
    config.add_argument("--filters", nargs="+",
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
        help="Show or hide the threshold by default. This will be ignored if no threshold is provided.")

    optional_settings = export.add_argument_group(
        title="OPTIONAL SETTINGS"
    )
    optional_settings.add_argument("--minify-json", action="store_true",
        help="Export JSON without indentation or line returns.")


    concat = subparsers.add_parser("concat", help="Concatenate multiple measurements JSONs into a single JSON file")
    concat.add_argument("--jsons", required=True, type=str, nargs="+", metavar="JSONs",
        help="Measurement JSON files to concatenate.")
    concat.add_argument("--default-collection", type=str,
        help="The key of the default collection to display. " +
             "If not provided, the first collection of the first JSON file will be displayed")
    concat.add_argument("--minify-json", action="store_true",
        help="Concatenate JSONs without indentation or line returns.")
    concat.add_argument("--output-json", required=True, metavar="JSON", type=str,
        help="Output JSON file")


def run(args):
    if args.subcommand == 'export':
        return export_measurements(vars(args))
    if args.subcommand == "concat":
        return concat_measurements(vars(args))
