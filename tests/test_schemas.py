import json
import jsonschema.validators
import pytest
from pathlib import Path

schemas = list(Path("augur/data/").glob("schema-*.json"))

@pytest.mark.parametrize("schema_path", schemas, ids = lambda schema_path: str(schema_path))
def test_schema_is_valid(schema_path):
    with schema_path.open("rb") as schema_fh:
        schema = json.load(schema_fh)

    Validator = jsonschema.validators.validator_for(schema)
    Validator.check_schema(schema)
