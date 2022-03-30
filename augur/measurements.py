"""
Create JSON files suitable for visualization within the measurements panel of Auspice.
"""
import os
import pandas as pd
import sys

from .utils import write_json
from .validate import measurements as validate_measurements_json, ValidateError


def column_exists(collection, column, column_purpose):
    """
    Checks the provided column exists in the provided collection.
    Prints an error message to stderr if the column does not exist.

    Parameters
    ----------
    collection: pandas.DataFrame
        Collection of measurements and metadata
    column: str
        Column to check exists in the collection
    column_purpose: str
        Purpose of provided column for detailed error message

    Returns
    -------
    bool
        True if column exists in collection
    """
    column_in_df = column in collection.columns
    if not column_in_df:
        print(
            f"ERROR: Provided {column_purpose} column «{column}» does not exist in collection TSV.",
            file=sys.stderr,
        )
    return column_in_df


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
            f"ERROR: collection TSV file ({collection}) does not exist",
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
        if not column_exists(collection_df, provided_column, required_column):
            checks_passed = False
        # Confirm the provided column does not overwrite an existing column
        if (required_column in collection_df.columns and
            provided_column != required_column):
            print(
                f"ERROR: Cannot use provided '{provided_column}' column as the {required_column} column " +
                f"because a '{required_column}' column already exists in collection TSV.",
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
        print(f"ERROR: Found a non-numeric measurement value: {e}", file=sys.stderr)
        return None

    return collection_df


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

    # Create collection output object with required keys
    collection_output = {
        'key': os.path.basename(args['collection']),
        'groupings': [{'key': col} for col in args['grouping_column'] if column_exists(collection_df, col, "grouping")],
        'x_axis_label': 'measurement values',
        'measurements': collection_df.to_dict(orient='records')
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


def register_arguments(parser):
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    export = subparsers.add_parser("export", help="Export a measurements JSON for a single collection")
    export.add_argument("--collection", required=True, metavar="TSV",
        help="Collection of measurements and metadata in a TSV file. " +
             "Keep in mind duplicate columns will be renamed as 'X', 'X.1', 'X.2'...'X.N'")
    export.add_argument("--strain-column", default="strain",
        help="Name of the column containing strain names. " +
             "Provided column will be renamed to `strain` so please make sure no other columns are named `strain`. " +
             "Strain names in this column should match the strain names in the corresponding Auspice dataset JSON.")
    export.add_argument("--value-column", default="value",
        help="Name of the column containing the numeric values to be plotted for the given collection. " +
             "Provided column will be renamed to `value` so please make sure no other columns are named `value`. ")
    export.add_argument("--grouping-column", required=True, nargs="+",
        help="Name of the column(s) that should be used as grouping(s) for measurements.")
    export.add_argument("--minify-json", action="store_true",
        help="Export JSON without indentation or line returns.")
    export.add_argument("--output-json", required=True, metavar="JSON", type=str,
        help="Output JSON file. " +
             "The file name must follow the Auspice sidecar file naming convention to be recognized as a sidecar file. " +
             "See Nextstrain data format docs for more details.")


def run(args):
    if args.subcommand == 'export':
        return export_measurements(vars(args))
