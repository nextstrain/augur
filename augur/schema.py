"""
Helpers for loading Augur JSON and YAML schemas.
"""
import json
from pathlib import Path
import jsonschema
import jsonschema.exceptions
from referencing import Registry
from ruamel.yaml import YAML, YAMLError

from augur.data import as_file
from augur.io.file import open_file


class InvalidSchemaError(Exception):
    pass


COMMAND_SCHEMAS = {
    "ancestral": "v1",
    "refine": "v1",
}


def load_augur_json_schema(schema: str):
    """
    Load a JSON schema from the augur included set of schemas
    (located in augur/data), specified by URL or local file name.
    """
    local_refs = _get_local_refs()
    if schema in local_refs:
        # Convert URL to local reference
        filename = local_refs[schema]
    else:
        # Assume everything else is already a local reference
        filename = schema

    try:
        with as_file(filename) as file:
            return load_json_schema_locally(file)
    except FileNotFoundError as err:
        raise FileNotFoundError(f"No local file for schema: {schema}") from err


def load_json_schema_locally(path):
    """
    Load a JSON schema offline using local files.

    If the schema references Augur schemas via URL, the Augur schemas are served
    by files bundled with the local Augur installation instead of what's served
    at the remote URL.

    There are two advantages to this approach:

    1. it uses the schema expected by the current Augur version (URL may serve
       schema for a newer version)
    2. it doesn't require an internet connection
    """
    return _load_json_schema(path, refs=_get_local_refs())


def _get_local_refs() -> dict[str, str]:
    """
    Return a mapping of Augur schema URL to local file name under augur/data.
    """
    local_refs = {
        'https://nextstrain.org/schemas/augur/annotations': "schema-annotations.json",
        'https://nextstrain.org/schemas/dataset/root-sequence': "schema-export-root-sequence.json",
        'https://nextstrain.org/schemas/auspice/config/v2': "schema-auspice-config-v2.json",
        'https://nextstrain.org/schemas/augur/subsample-config/v1': "schema-subsample-config.json",
        'https://nextstrain.org/schemas/augur/subsample-config-unaligned/v1': "schema-subsample-config-unaligned.json",
    }

    for command, version in COMMAND_SCHEMAS.items():
        command_name = command.replace(' ', '-')
        local_refs[f"https://nextstrain.org/schemas/augur/{command_name}-config/{version}"] = f"schema-{command_name}-config.json"

    return local_refs


def _load_json_schema(path, refs=None):
    """
    Load a schema from the path.

    Parameters
    ----------
    path
        Path to a JSON or YAML schema file.
    refs : dict
        Mapping of schema URI strings to local file paths (relative to
        augur/data) used to resolve schema references offline.
    """
    is_yaml = Path(path).suffix in ('.yaml', '.yml')
    try:
        with open_file(path, "r") as fh:
            if is_yaml:
                schema = YAML(typ="safe").load(fh)
            else:
                schema = json.load(fh)
    except (json.JSONDecodeError, YAMLError) as err:
        raise InvalidSchemaError(f"Schema {path} is not a valid {'YAML' if is_yaml else 'JSON'} file. Error: {err}")
    # check loaded schema is itself valid -- see http://python-jsonschema.readthedocs.io/en/latest/errors/
    Validator = jsonschema.validators.validator_for(schema)
    try:
        Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise InvalidSchemaError(f"Schema {path} is not a valid JSON Schema ({Validator.META_SCHEMA['$schema']}). Error: {err}")

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
