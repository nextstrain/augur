import pytest
import random

from augur.validate import (
    validate_collection_config_fields,
    validate_collection_display_defaults,
    validate_measurements_config,
    load_json_schema,
    validate_json,
    ValidateError
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


@pytest.fixture
def genome_annotation_schema():
    return load_json_schema("schema-annotations.json")

class TestValidateGenomeAnnotations():
    def test_negative_strand_nuc(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 200, "strand": "-"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_nuc_not_starting_at_one(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 100, "end": 200, "strand": "+"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_missing_nuc(self, capsys, genome_annotation_schema):
        d = {"cds": {"start": 100, "end": 200, "strand": "+"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_missing_properties(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 100}, "cds": {"start": 20, "strand": "+"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_not_stranded_cds(self, capsys, genome_annotation_schema):
        # Strand . is for features that are not stranded (as per GFF spec), and thus they're not CDSs
        d = {"nuc": {"start": 1, "end": 100}, "cds": {"start": 18, "end": 20, "strand": "."}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_negative_coordinates(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 100}, "cds": {"start": -2, "end": 10, "strand": "+"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_valid_genome(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 100}, "cds": {"start": 20,  "end": 28, "strand": "+"}}
        validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_valid_segmented_genome(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 100},
             "cds": {"segments": [{"start": 20,  "end": 28}], "strand": "+"}}
        validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_invalid_segmented_genome(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 100},
             "cds": {"segments": [{"start": 20,  "end": 28}, {"start": 27}], "strand": "+"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing

    def test_string_coordinates(self, capsys, genome_annotation_schema):
        d = {"nuc": {"start": 1, "end": 100},
             "cds": {"segments": [{"start": 20,  "end": 28}, {"start": "27", "end": "29"}], "strand": "+"}}
        with pytest.raises(ValidateError):
            validate_json(d, genome_annotation_schema, "<test-json>")
        capsys.readouterr() # suppress validation error printing