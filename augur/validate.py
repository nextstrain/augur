"""
Validate files related to augur consumption or export.
"""

import sys
import functools
import json
from .validate_export import (
    verifyMainJSONIsInternallyConsistent,
    verifyMetaAndOrTreeJSONsAreInternallyConsistent,
)

from augur.validation.json_validator import JsonValidator


class JsonValidationError(Exception):
    def __init__(self, errors):
        self.errors = errors


@functools.lru_cache()
def load_json(path):
    with open(path, "rb") as fh:
        return json.load(fh)


def validate_json_file(json_path, schema_name):
    print(f"Validating schema of {json_path} against {schema_name}...", end="")

    json_obj = load_json(json_path)
    json_validator = JsonValidator(json_obj, schema_name)

    if json_validator.is_valid:
        print("Passed")
        return

    for error in json_validator.errors:
        print(error, file=sys.stderr)

    raise JsonValidationError(json_validator.errors)


def auspice_config_v2(config_json, **_kwargs):
    validate_json_file(config_json, "schema-auspice-config-v2.json")


def export_v2(main_json, **_kwargs):
    if any(
        main_json.endswith(name)
        for name in {"frequencies.json", "entropy.json", "sequences.json"}
    ):
        raise JsonValidationError(
            "This validation subfunction is for the main `augur export v2` JSON only."
        )

    validate_json_file(main_json, "schema-export-v2.json")

    if not verifyMainJSONIsInternallyConsistent(
        load_json(main_json), JsonValidationError
    ):
        print(f"{main_json} is not internally-consistent.")


def export_v1(meta_json, tree_json, **_kwargs):
    if not meta_json.endswith("_meta.json"):
        raise JsonValidationError(
            f"The metadata JSON pathname {meta_json} must end with '_meta.json'."
        )

    if not tree_json.endswith("_tree.json"):
        raise JsonValidationError(
            f"The metadata JSON pathname {tree_json} must end with '_tree.json'."
        )

    validate_json_file(meta_json, "schema-export-v1-meta.json")
    validate_json_file(tree_json, "schema-export-v1-tree.json")

    if not verifyMetaAndOrTreeJSONsAreInternallyConsistent(
        load_json(meta_json), load_json(tree_json), JsonValidationError
    ):
        print(f"{meta_json} and {tree_json} are not internally-consistent.")


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
        globals()[args.subcommand.replace('-', '_')](**vars(args))
    except Exception as e:
        print(f"FATAL ERROR: {e}")
        sys.exit(2)
