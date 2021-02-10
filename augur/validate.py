"""
Validate files related to augur consumption or export.
"""

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
    json_obj = load_json(json_path)
    json_validator = JsonValidator(json_obj, schema_name)

    if not json_validator.is_valid:
        raise JsonValidationError(json_validator.errors)


def auspice_config_v2(fname):
    validate_json_file(fname, "schema-auspice-config-v2.json")


def export_v2(fname):
    if any(
        fname.endswith(name)
        for name in {"frequencies.json", "entropy.json", "sequences.json"}
    ):
        raise JsonValidationError(
            "This validation subfunction is for the main `augur export v2` JSON only."
        )

    validate_json_file(fname, "schema-export-v2.json")

    if not verifyMainJSONIsInternallyConsistent(load_json(fname), JsonValidationError):
        raise JsonValidationError(f"{fname} is not internally-consistent.")


def export_v1(meta_fname, tree_fname):
    if not meta_fname.endswith("_meta.json"):
        raise JsonValidationError(
            f"The metadata JSON pathname {meta_fname} must end with '_meta.json'."
        )

    if not tree_fname.endswith("_tree.json"):
        raise JsonValidationError(
            f"The metadata JSON pathname {tree_fname} must end with '_tree.json'."
        )

    validate_json_file(meta_fname, "schema-export-v1-meta.json")
    validate_json_file(tree_fname, "schema-export-v1-tree.json")

    if not verifyMetaAndOrTreeJSONsAreInternallyConsistent(
        load_json(meta_fname), load_json(tree_fname), JsonValidationError
    ):
        raise JsonValidationError(
            f"{meta_fname} and {tree_fname} are not internally-consistent."
        )
