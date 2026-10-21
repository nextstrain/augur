import pytest
import json
from ruamel.yaml import YAML

from augur.schema import (
    load_augur_json_schema,
    load_json_schema_locally,
    _load_json_schema,
)


def test_load_json_schema_locally(tmp_path):
    schema_file = tmp_path / "test-schema.json"
    raw_schema = {
        "type": "object",
        "properties": {
            "subsample": {
                "type": "object",
                "patternProperties": {
                    "^.*$": {
                        "$ref": "https://nextstrain.org/schemas/augur/subsample-config/v1"
                    }
                }
            }
        }
    }
    with open(schema_file, "w") as f:
        json.dump(raw_schema, f)

    validator = load_json_schema_locally(schema_file)
    valid_data = {
        "subsample": {
            "global": {
                "samples": {
                    "focal": {
                        "group_by": ["country"]
                    }
                }
            }
        }
    }
    validate_json(valid_data, validator, "")


def test_load_augur_json_schema():
    validator = load_augur_json_schema("schema-annotations.json")
    assert validator.schema["$id"] == "https://nextstrain.org/schemas/augur/annotations"

    validator = load_augur_json_schema("https://nextstrain.org/schemas/augur/tree-config/v1")
    assert validator.schema["$id"] == "https://nextstrain.org/schemas/augur/tree-config/v1"

    with pytest.raises(FileNotFoundError, match="No local file for schema"):
        load_augur_json_schema("https://nextstrain.org/schemas/augur/invalid-config/v1")


@pytest.mark.parametrize("ext", ["json", "yaml", "yml"])
def test_load_json_schema(tmp_path, ext):
    schema_file = tmp_path / f"test-schema.{ext}"
    raw_schema = {
        "type": "object",
        "properties": {
            "name": {
                "type": "string"
            }
        }
    }
    with open(schema_file, "w") as f:
        if ext == "json":
            json.dump(raw_schema, f)
        else:
            YAML(typ="safe").dump(raw_schema, f)

    validator = _load_json_schema(schema_file)
    valid_data = {"name": "test"}
    validate_json(valid_data, validator, "")
