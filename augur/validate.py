"""
Validate files related to augur consumption or export.
"""

import sys
from collections import defaultdict
import json
import jsonschema
import jsonschema.exceptions
import re
from itertools import groupby
from referencing import Registry
from textwrap import indent
from typing import Iterable, Union
from augur.data import as_file
from augur.io.file import open_file
from augur.io.print import print_err
from augur.io.json import shorten_as_json
from .validate_export import verifyMainJSONIsInternallyConsistent, verifyMetaAndOrTreeJSONsAreInternallyConsistent
from .types import ValidationMode

def fatal(message):
    print_err("FATAL ERROR: {}".format(message))
    sys.exit(2)

class ValidateError(Exception):
    pass

def validation_failure(mode: ValidationMode):
    if mode is ValidationMode.ERROR:
        sys.exit(2)
    elif mode is ValidationMode.WARN:
        print_err(f"Continuing due to --validation-mode={mode.value} even though there were validation errors.")
    elif mode is ValidationMode.SKIP:
        # Shouldn't be doing validation under skip, but if we're called anyway just do nothing.
        return
    else:
        raise ValueError(f"unknown validation mode: {mode!r}")

def load_json_schema(path, refs=None):
    '''
    Load a JSON schema from the augur included set of schemas
    (located in augur/data)
    '''
    try:
        with as_file(path) as file, open_file(file, "r") as fh:
            schema = json.load(fh)
    except json.JSONDecodeError as err:
        raise ValidateError("Schema {} is not a valid JSON file. Error: {}".format(path, err))
    # check loaded schema is itself valid -- see http://python-jsonschema.readthedocs.io/en/latest/errors/
    Validator = jsonschema.validators.validator_for(schema)
    try:
        Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise ValidateError(f"Schema {path} is not a valid JSON Schema ({Validator.META_SCHEMA['$schema']}). Error: {err}")

    if refs:
        # Make the validator aware of additional schemas
        schema_store = dict()
        for k, v in refs.items():
            with as_file(v) as file, open_file(file, "r") as fh:
                schema_store[k] = json.load(fh)

        # Create a dummy retrieval function to handle URIs not present in
        # schema_store. This often indicates a typo (the $ref doesn't match the
        # key of the schema_store) or we forgot to add a local mapping for a new
        # $ref.
        def retrieve(uri):
            # Take advantage of the fact that BaseException is not handled by
            # Registry.get_or_retrieve. This means the custom error message is
            # printed instead of the less helpful default:
            #    jsonschema.exceptions._WrappedReferencingError: Unresolvable: https://…
            raise BaseException(f"The schema used for validation could not resolve a local file for {uri!r}. " +
                            "Please check the schema used and update the appropriate schema_store as needed." )

        registry = Registry(retrieve=retrieve).with_contents(schema_store.items())
        schema_validator = Validator(schema, registry=registry)
    else:
        schema_validator = Validator(schema)

    return schema_validator


def load_json(path):
    with open_file(path, 'rb') as fh:
        try:
            jsonToValidate = json.load(fh)
        except json.JSONDecodeError:
            raise ValidateError("Supplied JSON to validate ({}) is not a valid JSON".format(path))
    return jsonToValidate

def validate_json(jsonToValidate, schema, filename):
    # See <https://python-jsonschema.readthedocs.io/en/v3.2.0/errors/> and
    # <https://python-jsonschema.readthedocs.io/en/v3.2.0/validate/> for the
    # jsonschema APIs we use here.
    print_err("Validating schema of {!r}...".format(filename))

    # Find all errors.  This is what schema.validate() uses internally, before
    # raising just one "best" error.  We want to report ~everything at once, so
    # the user isn't stuck playing whack-a-mole.
    errors = list(schema.iter_errors(jsonToValidate))

    def custom_message(error):
        """Convert technical JSON schema errors to human-readable messages."""
        validator = error.validator
        short_value = shorten_as_json(error.instance, 50, "…")

        if validator == "oneOf":
            return f"{short_value} did not match one of the acceptable options below."
        if validator == "anyOf":
            return f"{short_value} did not match any of the acceptable options below."

        if validator == "required":
            missing_property = str(error.args[0]).split("'")[1]
            return f"Missing required property '{missing_property}'"

        if validator == "additionalProperties":
            additional_property = list(error.instance.keys() - error.schema.get('properties', {}).keys())[0]
            return f"Unexpected property '{additional_property}'"

        if validator == "type":
            expected_type = error.validator_value
            actual_type = type(error.instance).__name__
            return f"Expected {expected_type} but found {actual_type} {short_value!r}"

        if validator == "const":
            expected_value = error.validator_value
            return f"Expected '{expected_value}' but found {error.instance!r}"

        if validator == "pattern":
            pattern = error.validator_value
            if pattern == "^[0-9X]{4}-[0-9X]{2}-[0-9X]{2}$":
                return f"Expected a date in the format YYYY-MM-DD but found {short_value}"
            else:
                return f"Expected a value with the pattern {pattern} but found {short_value}"

        # Fallback to an error message from jsonschema.
        if error.args:
            return error.args[0]

        return "Unknown validation error"

    def print_errors(errors, level=1):
        prefix = "  "

        for error in sorted(errors, key=jsonschema.exceptions.relevance):
            path = elide_path(error.absolute_path)

            print_err(indent(f"{path or 'top level'} failed: {custom_message(error)}", prefix*level))

            # Report sub-errors, as they're often closer to what needs fixing.
            #
            # We special case the oneOf and anyOf validators to partition
            # errors by the validation arm, so it's clear why every arm of the
            # oneOf/anyOf validator failed.  This is particularly useful for
            # our top-level "tree" property, which uses oneOf.
            #
            # In theory this special handling could also apply to allOf, not,
            # if/then/else, items, etc. validators which take subschemas, but
            # jsonschema 3.2.0 doesn't include context for those.  We currently
            # only use oneOf any how.
            #   -trs, 23 Jan 2023
            if error.context:
                if error.validator in {"oneOf", "anyOf"}:
                    validator_value_idx = lambda e: e.schema_path[0]
                    for idx, ctx in grouped(error.context, key=validator_value_idx):
                        validator_subvalue = shorten_as_json(error.validator_value[idx], 100, "…")
                        print_err(indent(f"Option {idx+1}: {validator_subvalue}", prefix*(level+1)))
                        print_errors(ctx, level+2)
                else:
                    print_errors(error.context, level+1)

    print_errors(errors)

    if errors:
        raise ValidateError("Validation of {!r} failed.".format(filename))

def elide_path(path: Iterable[Union[str, int]]) -> str:
    """
    Elide recursive ``.children[N]`` selectors from *path* after formatting
    with :py:func:`format_path`.

    Our tree JSONs are highly nested and the exact path to a node in the tree
    is often more noisy than useful in error messages.

    >>> format_path(("tree", "children", 1, "children", 0, "children", 1, "children", 0, "node_attrs", "x"))
    '.tree.children[1].children[0].children[1].children[0].node_attrs.x'
    >>> elide_path(("tree", "children", 1, "children", 0, "children", 1, "children", 0, "node_attrs", "x"))
    '.tree.children[…].node_attrs.x'

    Elision is only used with more than a single level of ``.children[N]``.

    >>> elide_path(("tree", "children", 1, "node_attrs", "x"))
    '.tree.children[1].node_attrs.x'
    """
    return re.sub(r'([.]children\[[0-9]+\]){2,}', r'.children[…]', format_path(path))

def format_path(path: Iterable[Union[str, int]]) -> str:
    """
    Format an iterable of *path* segments, which index into a JSON document,
    into a more human-readable string.

    Intended for folks who aren't necessarily programmers, so this doesn't try
    to construct a valid JS property accessor chain or valid JSON Path
    selector.

    >>> format_path(("a", "b", "c"))
    '.a.b.c'
    >>> format_path(("l", "m-n-o", "p"))
    ".l.'m-n-o'.p"
    >>> format_path(("x", "y", 42, "z"))
    '.x.y[42].z'
    """
    def valid_identifier(x) -> bool:
        return isinstance(x, str) and re.search(r'^[a-zA-Z$_][a-zA-Z0-9_$]*$', x) is not None

    def fmt(x) -> str:
        return (f"[{x}]"  if isinstance(x, int)  else
                f".{x}"   if valid_identifier(x) else
                f".{x!r}")

    return "".join(map(fmt, path))

def grouped(iterable, key):
    """
    Version of :py:func:`itertools.groupby` which doesn't require the caller to
    remember to sort first.
    """
    return groupby(sorted(iterable, key=key), key=key)


validate = validate_json  # TODO update uses and drop this alias


def auspice_config_v2(config_json: Union[str,dict], **kwargs):
    schema = load_json_schema("schema-auspice-config-v2.json")
    if isinstance(config_json, dict):
        config = config_json
        filename = "merged config"
    else:
        config = load_json(config_json)
        filename = config_json

    validate(config, schema, filename)

def export_v2(main_json, **kwargs):
    # The main_schema uses references to other schemas, and the suggested use is
    # to define these refs as valid URLs. Augur itself should not access schemas
    # over the wire so we provide a mapping between URLs and filepaths here. The
    # filepath is specified relative to ./augur/data (where all the schemas
    # live).
    refs = {
        'https://nextstrain.org/schemas/augur/annotations': "schema-annotations.json",
        'https://nextstrain.org/schemas/dataset/root-sequence': "schema-export-root-sequence.json",
        'https://nextstrain.org/schemas/auspice/config/v2': "schema-auspice-config-v2.json",
    }
    main_schema = load_json_schema("schema-export-v2.json", refs)

    if main_json.endswith("frequencies.json") or main_json.endswith("entropy.json") or main_json.endswith("sequences.json"):
        raise ValidateError("This validation subfunction is for the main `augur export v2` JSON only.")

    main = load_json(main_json)
    validate(main, main_schema, main_json)

    if verifyMainJSONIsInternallyConsistent(main, ValidateError):
        print_err("Validation of {!r} succeeded.".format(main_json))
    else:
        print_err("Validation of {!r} succeeded, but there were warnings you may want to resolve.".format(main_json))


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
        print_err("Validation of {!r} and {!r} succeeded.".format(meta_json, tree_json))
    else:
        print_err("Validation of {!r} and {!r} succeeded, but there were warnings you may want to resolve.".format(meta_json, tree_json))


def get_unique_keys(list_of_dicts):
    """
    Returns a set of unique keys from a list of dicts

    Examples
    --------
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
            print_err(
                f"ERROR: Collection {include_index}includes {config_field} that",
                f"do not exist as fields in measurements: {invalid_fields}."
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
        print_err(
            f"ERROR: Collection {include_index}has a default group-by field",
            f"'{default_grouping}' that is not included in the groupings' fields."
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
            print_err(
                f"ERROR: Collections at indexes {collection_indexes} share the same collection key '{collection_key}'."
            )

    # Check the default collection value matches a collection's key value
    default_collection = measurements.get('default_collection')
    if default_collection and default_collection not in collection_keys.keys():
        valid_measurements_config = False
        print_err(
            f"ERROR: The default collection key {default_collection!r} does not match any of the collections' keys."
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
