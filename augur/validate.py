"""
Validate files related to augur consumption or export.
"""

import sys
import os
from collections import defaultdict
import json
import jsonschema
from pkg_resources import resource_string
from .validate_export import verifyMainJSONIsInternallyConsistent, verifyMetaAndOrTreeJSONsAreInternallyConsistent

class ValidationWarnings:
    def __init__(self):
        self.seen = defaultdict(set)
    def add(self, warningType, message):
        self.seen[warningType].add(message)
    def show(self):
        print("WARNINGS")
        print(self.seen)

class ValidationErrors(ValidationWarnings):
    def show(self):
        print("ERRORS")
        print(self.seen)
        sys.exit(2)

def fatal(message):
    print("FATAL ERROR: {}".format(message))
    sys.exit(2)

class ValidateError(Exception):
    pass


def load_json_schema(path):
    '''
    Load a JSON schema from the augur included set of schemas
    (located in augur/data)
    '''
    try:
        schema = json.loads(resource_string(__package__, os.path.join("data", path)))
    except json.JSONDecodeError as err:
        raise ValidateError("Schema {} is not a valid JSON file. Error: {}".format(path, err))
    # check loaded schema is itself valid -- see http://python-jsonschema.readthedocs.io/en/latest/errors/
    Validator = jsonschema.validators.validator_for(schema)
    try:
        Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise ValidateError(f"Schema {path} is not a valid JSON Schema ({Validator.META_SCHEMA['$schema']}). Error: {err}")
    return Validator(schema)

def load_json(path):
    with open(path, 'rb') as fh:
        try:
            jsonToValidate = json.load(fh)
        except json.JSONDecodeError:
            raise ValidateError("Supplied JSON to validate ({}) is not a valid JSON".format(path))
    return jsonToValidate

def validate_json(jsonToValidate, schema, filename):
    print("Validating schema of {!r}...".format(filename))
    try:
        schema.validate(jsonToValidate) # https://python-jsonschema.readthedocs.io/en/latest/validate/
    except jsonschema.exceptions.ValidationError:
        for error in sorted(schema.iter_errors(jsonToValidate), key=str):
            trace = list(error.schema_path) # this may be really long for nested tree structures
            if len(trace) > 6:
                trace = ["..."] + trace[-5:]
            trace = [str(x) for x in trace]
            print("\tERROR: {}. Trace: {}".format(error.message, " - ".join(trace)), file=sys.stderr)
        raise ValidateError("Validation of {!r} failed.".format(filename))


validate = validate_json  # TODO update uses and drop this alias


def auspice_config_v2(config_json, **kwargs):
    schema = load_json_schema("schema-auspice-config-v2.json")
    config = load_json(config_json)
    validate(config, schema, config_json)

def export_v2(main_json, **kwargs):
    # validationWarnings = ValidationWarnings()
    # validationErrors = ValidationErrors()
    main_schema = load_json_schema("schema-export-v2.json")

    if main_json.endswith("frequencies.json") or main_json.endswith("entropy.json") or main_json.endswith("sequences.json"):
        raise ValidateError("This validation subfunction is for the main `augur export v2` JSON only.")

    main = load_json(main_json)
    validate(main, main_schema, main_json)

    if verifyMainJSONIsInternallyConsistent(main, ValidateError):
        print("Validation of {!r} succeeded.".format(main_json))
    else:
        print("Validation of {!r} succeeded, but there were warnings you may want to resolve.".format(main_json))


def export_v1(meta_json, tree_json, **kwargs):
    meta_schema = load_json_schema("schema-export-v1-meta.json")
    tree_schema = load_json_schema("schema-export-v1-tree.json")

    if not meta_json.endswith("_meta.json"):
        raise ValidateError("The metadata JSON pathname {} must end with '_meta.json'.".format(meta_json))

    if not tree_json.endswith("_tree.json"):
        raise ValidateError("The metadata JSON pathname {} must end with '_tree.json'.".format(tree_json))

    meta = load_json(meta_json)
    tree = load_json(tree_json)

    validate(meta, meta_schema, meta_json)
    validate(tree, tree_schema, tree_json)

    if verifyMetaAndOrTreeJSONsAreInternallyConsistent(meta, tree, ValidateError):
        print("Validation of {!r} and {!r} succeeded.".format(meta_json, tree_json))
    else:
        print("Validation of {!r} and {!r} succeeded, but there were warnings you may want to resolve.".format(meta_json, tree_json))


def get_unique_keys(list_of_dicts):
    """
    Returns a set of unique keys from a list of dicts

    >>> list_of_dicts = [{"key1": "val1", "key2": "val2"}, {"key1": "val1", "key3": "val3"}]
    >>> sorted(get_unique_keys(list_of_dicts))
    ['key1', 'key2', 'key3']
    """
    return set().union(*(single_dict.keys() for single_dict in list_of_dicts))


def validate_collection_config_fields(collection, index=None):
    """
    Validates a single collection's config field keys provided in fields,
    groupings, and filters are valid fields that exist in measurements' fields.

    Prints any validation errors to stderr.

    Parameters
    ----------
    collection: dict
        A single collection to validate. Assumes that the collection has already passed the schema validation.
    index: int, optional
        the index of the collection within a list of collections in a measurements JSON.
        Used to print more detailed error messages.

    Returns
    -------
    bool
        True if collection's config is valid
    """
    valid_collection_config_fields = True
    nested_config_fields = ['fields', 'groupings']
    flat_config_fields = ['filters']
    # Create set of all measurements' fields for verifying field configs
    all_measurement_fields = get_unique_keys(collection['measurements'])

    for config_field in (nested_config_fields + flat_config_fields):
        invalid_fields = set()
        for config_value in collection.get(config_field, []):
            # config value can be a field name string (i.e. flat_config_fields)
            # or a dict with the field name in 'key' (i.e. nested_config_fields)
            field_name = config_value['key'] if config_field in nested_config_fields else config_value
            if field_name not in all_measurement_fields:
                invalid_fields.add(field_name)

        if invalid_fields:
            valid_collection_config_fields = False
            include_index = f"(at index {index}) " if index is not None else ""
            print(
                f"ERROR: Collection {include_index}includes {config_field} that",
                f"do not exist as fields in measurements: {invalid_fields}.",
                file=sys.stderr
            )

    return valid_collection_config_fields


def validate_collection_display_defaults(collection, index=None):
    """
    Validates a single collection's display defaults. If a default group-by
    field is provided, the field must be included in groupings.

    Prints validation errors to stderr.

    Parameters
    ----------
    collection: dict
        A single collection to validate. Assumest htat the collection has already passed the schema validation.
    index: int, optional
        The index of the collection within a list of collections in a measurements JSON.
        Used to print more detailed error messages.

    Returns
    -------
    bool
        True if collection's display defaults are valid
    """
    valid_display_defaults = True

    grouping_fields = {grouping['key'] for grouping in collection['groupings']}
    default_grouping = collection.get('display_defaults', {}).get('group_by')

    if default_grouping and default_grouping not in grouping_fields:
        valid_display_defaults = False
        include_index = f"(at index {index}) " if index is not None else ""
        print(
            f"ERROR: Collection {include_index}has a default group-by field",
            f"'{default_grouping}' that is not included in the groupings' fields.",
            file=sys.stderr
        )

    return valid_display_defaults


def validate_measurements_config(measurements):
    """
    Validate measurements' config values meet expectations described in the
    measurements JSON schema descriptions that cannot be verified via
    `validate_json`:
    1. Individual collections have valid config values
    2. All collections have unique keys
    3. If a default collection is provided, it matches one of the collections

    Prints any validation errors to stderr.

    Parameters
    ----------
    measurements: dict
        Loaded measurements JSON to validate. Assumes the measurements JSON has already passed the schema validation.

    Returns
    -------
    bool
        True if measurements' config is valid
    """
    valid_measurements_config = True
    collection_keys = defaultdict(list)

    # First check configs for individual collections
    for index, collection in enumerate(measurements['collections']):
        # Save the collection key and index of collection to verify unique keys later
        collection_keys[collection['key']].append(index)

        if not all([
            validate_collection_config_fields(collection, index),
            validate_collection_display_defaults(collection, index)
        ]):
            valid_measurements_config = False

    # Check collections have unique keys
    for collection_key, collection_indexes in collection_keys.items():
        if len(collection_indexes) > 1:
            valid_measurements_config = False
            print(
                f"ERROR: Collections at indexes {collection_indexes} share the same collection key '{collection_key}'.",
                file=sys.stderr
            )

    # Check the default collection value matches a collection's key value
    default_collection = measurements.get('default_collection')
    if default_collection and default_collection not in collection_keys.keys():
        valid_measurements_config = False
        print(
            f"ERROR: The default collection key {default_collection!r} does not match any of the collections' keys.",
            file=sys.stderr
        )

    return valid_measurements_config


def measurements(measurements_json, **kwargs):
    schema = load_json_schema("schema-measurements.json")
    measurements = load_json(measurements_json)
    validate_json(measurements, schema, measurements_json)
    if not validate_measurements_config(measurements):
        raise ValidateError("Validation of the measurements' config values failed.")
    return measurements


def measurements_collection_config(collection_config_json, **kwargs):
    schema = load_json_schema("schema-measurements-collection-config.json")
    collection_config = load_json(collection_config_json)
    validate_json(collection_config, schema, collection_config_json)
    if not validate_collection_display_defaults(collection_config):
        raise ValidateError("Validation of the collection config display defaults failed.")
    return collection_config


def register_parser(parent_subparsers):
    # Not using utils.first_line for help here because it results in a circular import
    parser = parent_subparsers.add_parser("validate", help=__doc__)
    subparsers = parser.add_subparsers(dest="subcommand", help="Which file(s) do you want to validate?")

    subparsers.add_parser("export-v2", help="validate JSON intended for auspice v2") \
        .add_argument('main_json', metavar='JSON', help="exported (main) v2 auspice JSON")

    export_v1 = subparsers.add_parser("export-v1", help="validate tree+meta JSONs intended for auspice v1")
    export_v1.add_argument('meta_json', metavar='META-JSON', help="exported (v1) meta JSON")
    export_v1.add_argument('tree_json', metavar='TREE-JSON', help="exported (v1) tree JSON")

    subparsers.add_parser("auspice-config-v2", help="validate auspice config intended for `augur export v2`") \
        .add_argument('config_json', metavar='JSON', help="auspice config JSON")

    subparsers.add_parser("measurements", help="validate measurements JSON intended for auspice measurements panel") \
        .add_argument("measurements_json", metavar="JSON", help="exported measurements JSON")

    subparsers.add_parser("measurements-collection-config", help="validate measurement collection config intended for `augur measurements export`") \
        .add_argument("collection_config_json", metavar="JSON", help="collection config JSON")
    return parser


def run(args):
    try:
        globals()[args.subcommand.replace('-','_')](**vars(args))
    except ValidateError as e:
        fatal(e)
