import pytest
import random

from augur.validate import (
    validate_collection_config_fields,
    validate_collection_display_defaults,
    validate_measurements_config
)


@pytest.fixture
def example_collection_measurements():
    return [
        {"strain": "strain_1", "value": 0, "valid_field_1": "value_1a", "valid_field_2": "value_2a", "valid_field_3": "value_3a"},
        {"strain": "strain_2", "value": 0, "valid_field_1": "value_1b", "valid_field_2": "value_2b", "valid_field_3": "value_3b"},
        {"strain": "strain_3", "value": 0, "valid_field_1": "value_1c", "valid_field_2": "value_2c", "valid_field_3": "value_3c"}
    ]

@pytest.fixture
def example_collection(example_collection_measurements):
    return {
        "key": "collection_0",
        "fields": [
            {"key": "valid_field_1"},
            {"key": "valid_field_2"},
            {"key": "valid_field_3"}
        ],
        "groupings": [
            {"key": "valid_field_1"},
            {"key": "valid_field_2"},
            {"key": "valid_field_3"}
        ],
        "filters": ["valid_field_1", "valid_field_2", "valid_field_3"],
        "display_defaults": {
            "group_by": "valid_field_1"
        },
        "measurements": example_collection_measurements
    }

@pytest.fixture
def example_measurements(example_collection):
    number_of_collections = 10
    return {
        "default_collection": f"collection_{random.randint(0, number_of_collections - 1)}",
        "collections": [{**example_collection, "key": f"collection_{x}"} for x in range(number_of_collections)]
    }

class TestValidateMeasurements():
    def test_validate_collection_config_fields_valid(self, example_collection):
        assert validate_collection_config_fields(example_collection)

    @pytest.mark.parametrize(
        "invalid_config",
        [
            {"fields": [{"key": "invalid_field"}]},
            {"groupings": [{"key": "invalid_field"}]},
            {"filters": ["invalid_field"]}
        ]
    )
    def test_validate_collection_config_fields_invalid(self, invalid_config, example_collection_measurements, capsys):
        collection = {**invalid_config, "measurements": example_collection_measurements}
        assert not validate_collection_config_fields(collection)
        assert capsys.readouterr().err == f"ERROR: Collection includes {next(iter(invalid_config))} that do not exist as fields in measurements: {{'invalid_field'}}.\n"

    def test_validate_collection_display_defaults_valid(self, example_collection):
        assert validate_collection_display_defaults(example_collection)

    def test_validate_collection_display_defaults_invalid(self, example_collection, capsys):
        collection = {**example_collection}
        collection["display_defaults"]["group_by"] = "invalid_field"
        assert not validate_collection_display_defaults(collection)
        assert capsys.readouterr().err == "ERROR: Collection has a default group-by field 'invalid_field' that is not included in the groupings' fields.\n"

    def test_validate_measurements_config_valid(self, example_measurements):
        assert validate_measurements_config(example_measurements)

    def test_validate_measurements_config_duplicate_collection_keys(self, example_collection, capsys):
        measurements = {
            "collections": [example_collection] * 2
        }
        assert not validate_measurements_config(measurements)
        assert capsys.readouterr().err == "ERROR: Collections at indexes [0, 1] share the same collection key 'collection_0'.\n"

    def test_validate_measurements_config_invalid_default_collection(self, example_measurements, capsys):
        measurements = {
            **example_measurements,
            "default_collection": "invalid_collection"
        }
        assert not validate_measurements_config(measurements)
        assert capsys.readouterr().err == "ERROR: The default collection key 'invalid_collection' does not match any of the collections' keys.\n"
