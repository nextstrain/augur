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
    try:
        jsonschema.Draft6Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise ValidateError("Schema {} is not a valid JSON file. Error: {}".format(path, err))
    return jsonschema.Draft6Validator(schema)

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


def register_arguments(parser):
    subparsers = parser.add_subparsers(dest="subcommand", help="Which file(s) do you want to validate?")

    subparsers.add_parser("export-v2", help="validate JSON intended for auspice v2") \
        .add_argument('main_json', metavar='JSON', help="exported (main) v2 auspice JSON")

    export_v1 = subparsers.add_parser("export-v1", help="validate tree+meta JSONs intended for auspice v1")
    export_v1.add_argument('meta_json', metavar='META-JSON', help="exported (v1) meta JSON")
    export_v1.add_argument('tree_json', metavar='TREE-JSON', help="exported (v1) tree JSON")

    subparsers.add_parser("auspice-config-v2", help="validate auspice config intended for `augur export v2`") \
        .add_argument('config_json', metavar='JSON', help="auspice config JSON")


def run(args):
    try:
        globals()[args.subcommand.replace('-','_')](**vars(args))
    except ValidateError as e:
        fatal(e)
