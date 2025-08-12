"""Tests for augur subsample

This file contains unit tests for the augur subsample command.
"""

import pytest
import yaml
from textwrap import dedent

from augur import subsample
from augur.errors import AugurError


# Test fixtures for common data and mocks
@pytest.fixture
def valid_config_dict():
    """A valid subsample configuration as a Python dictionary."""
    return {
        "samples": {
            "focal": {
                "min_length": 8000,
                "group_by": ["region"],
                "sequences_per_group": 10,
                "probabilistic_sampling": True
            },
            "contextual": {
                "exclude_where": ["region=North America"],
                "min_date": "2020-01-01",
                "max_date": "2021-12-31",
                "subsample_max_sequences": 100
            }
        }
    }


@pytest.fixture
def valid_config_yaml(valid_config_dict):
    """A valid subsample configuration as a YAML string."""
    return yaml.dump(valid_config_dict)


@pytest.fixture
def config_with_root():
    """Configuration nested under a top level key."""
    return {
        "subsample": {
            "samples": {
                "test_sample": {
                    "min_length": 8000,
                }
            }
        },
        "other_data": "ignored"
    }


class TestParseConfig:
    """Tests for the _parse_config function."""

    def _create_config_file(self, tmp_path, content, filename="config.yaml"):
        """Helper to create a temporary config file."""
        config_file = tmp_path / filename
        config_file.write_text(content)
        return str(config_file)

    def test_parse_valid_yaml_config(self, valid_config_dict, valid_config_yaml, tmp_path):
        """Test parsing a valid YAML configuration file."""
        config_path = self._create_config_file(tmp_path, valid_config_yaml)
        result = subsample._parse_config(config_path)
        assert result == valid_config_dict

    def test_timestamps_as_strings(self, tmp_path):
        """Test that timestamp-like values are treated as strings, not parsed as dates."""
        yaml_with_timestamp = dedent("""\
            samples:
              test:
                min_date: 2020-03-01
                max_date: 2021-12-31
        """)
        config_path = self._create_config_file(tmp_path, yaml_with_timestamp)
        result = subsample._parse_config(config_path)
        assert isinstance(result["samples"]["test"]["min_date"], str)
        assert isinstance(result["samples"]["test"]["max_date"], str)
        assert result["samples"]["test"]["min_date"] == "2020-03-01"
        assert result["samples"]["test"]["max_date"] == "2021-12-31"

    def test_config_root_happy_path(self, config_with_root, tmp_path):
        """Test config-root extraction when the key exists."""
        yaml_data = yaml.dump(config_with_root)
        config_path = self._create_config_file(tmp_path, yaml_data)
        result = subsample._parse_config(config_path, "subsample")
        assert result == config_with_root["subsample"]

    def test_config_root_missing_key(self, config_with_root, tmp_path):
        """Test config-root extraction when the key doesn't exist."""
        yaml_data = yaml.dump(config_with_root)
        config_path = self._create_config_file(tmp_path, yaml_data)
        with pytest.raises(AugurError) as exc_info:
            subsample._parse_config(config_path, "missing_key")
        assert f"Config root key 'missing_key' not found in '{config_path}'" in str(exc_info.value)

    def test_yaml_syntax_error(self, tmp_path):
        """Test handling of YAML syntax errors."""
        invalid_yaml = "samples:\n  test:\n    invalid: [\n"  # Unclosed bracket
        config_path = self._create_config_file(tmp_path, invalid_yaml)
        with pytest.raises(AugurError) as exc_info:
            subsample._parse_config(config_path)
        assert f"Error parsing subsampling scheme '{config_path}'" in str(exc_info.value)

    def test_schema_validation_success(self, valid_config_yaml, tmp_path):
        """Test successful schema validation."""
        config_path = self._create_config_file(tmp_path, valid_config_yaml)
        result = subsample._parse_config(config_path)
        assert "samples" in result

    def test_schema_validation_failure(self, tmp_path):
        """Test handling of schema validation errors with invalid config."""
        invalid_yaml = dedent("""\
            # Missing required 'samples' key
            invalid_key: "test"
        """)
        config_path = self._create_config_file(tmp_path, invalid_yaml)
        with pytest.raises(AugurError) as exc_info:
            subsample._parse_config(config_path)
        assert "Config validation failed:" in str(exc_info.value)


class TestAddToArgs:
    """Tests for the _add_to_args helper function."""

    def test_scalar_value(self):
        """Test adding scalar values to arguments dict."""
        args = {}
        subsample._add_to_args(args, "--test-flag", "test_value")
        assert args == {"--test-flag": "test_value"}

    def test_list_value(self):
        """Test adding list values to arguments dict."""
        args = {}
        subsample._add_to_args(args, "--test-flag", ["val1", "val2", "val3"])
        assert args == {"--test-flag": ["val1", "val2", "val3"]}

    def test_tuple_value(self):
        """Test adding tuple values to arguments dict."""
        args = {}
        subsample._add_to_args(args, "--test-flag", ("val1", "val2"))
        assert args == {"--test-flag": ("val1", "val2")}

    def test_boolean_true_with_two_flags(self):
        """Test boolean true with tuple of (true_flag, false_flag)."""
        args = {}
        subsample._add_to_args(args, ("--enable", "--disable"), True)
        assert args == {"--enable": None}

    def test_boolean_false_with_false_flag(self):
        """Test boolean false with tuple containing false flag."""
        args = {}
        subsample._add_to_args(args, ("--enable", "--disable"), False)
        assert args == {"--disable": None}

    def test_boolean_false_with_no_false_flag(self):
        """Test boolean false with tuple containing None as false flag."""
        args = {}
        subsample._add_to_args(args, ("--enable", None), False)
        assert args == {}  # Nothing should be added

    def test_multiple_calls_preserve_flags(self):
        """Test that multiple calls add all flags to the dict."""
        args = {}
        subsample._add_to_args(args, "--first", "value1")
        subsample._add_to_args(args, "--second", ["val2", "val3"])
        subsample._add_to_args(args, ("--flag", None), True)
        expected = {
            "--first": "value1",
            "--second": ["val2", "val3"],
            "--flag": None
        }
        assert args == expected


class TestSampleConstructFilterArgs:
    """Tests for Sample._construct_filter_args method."""

    def test_always_includes_required_flags(self):
        """Test that required flags are always included."""
        sample = subsample.Sample("test_sample", {}, {})

        # Check required flags are present
        assert "--skip-checks" in sample.filter_args
        assert "--nthreads" in sample.filter_args
        assert "--output-strains" in sample.filter_args

        # Check values
        assert sample.filter_args["--skip-checks"] is None  # Boolean flag
        assert sample.filter_args["--nthreads"] == 1
        assert sample.filter_args["--output-strains"] == sample.output_strains

    def test_extends_global_args(self):
        """Test that global filter args are extended."""
        global_args = {"--metadata": "test.tsv", "--sequences": "test.fasta"}
        sample = subsample.Sample("test_sample", {}, global_args)

        # All global args should be present
        for flag, value in global_args.items():
            assert flag in sample.filter_args
            assert sample.filter_args[flag] == value

    def test_maps_yaml_keys_to_filter_flags(self):
        """Test mapping of YAML config keys to filter flags."""
        config = {
            "exclude": ["file1.txt", "file2.txt"],
            "include_where": ["region=Europe"],
            "min_date": "2020-01-01",
            "group_by": ["region", "country"],
            "sequences_per_group": 5,
            "subsample_max_sequences": 100
        }

        sample = subsample.Sample("test_sample", config, {})

        # Check flag mappings and values
        assert "--exclude" in sample.filter_args
        assert sample.filter_args["--exclude"] == ["file1.txt", "file2.txt"]
        assert "--include-where" in sample.filter_args
        assert sample.filter_args["--include-where"] == ["region=Europe"]
        assert "--min-date" in sample.filter_args
        assert sample.filter_args["--min-date"] == "2020-01-01"
        assert "--group-by" in sample.filter_args
        assert sample.filter_args["--group-by"] == ["region", "country"]
        assert "--sequences-per-group" in sample.filter_args
        assert sample.filter_args["--sequences-per-group"] == 5
        assert "--subsample-max-sequences" in sample.filter_args
        assert sample.filter_args["--subsample-max-sequences"] == 100

    def test_boolean_mapping_true(self):
        """Test boolean true mapping for probabilistic_sampling."""
        config = {"probabilistic_sampling": True}
        sample = subsample.Sample("test_sample", config, {})

        assert "--probabilistic-sampling" in sample.filter_args
        assert sample.filter_args["--probabilistic-sampling"] is None
        assert "--no-probabilistic-sampling" not in sample.filter_args

    def test_boolean_mapping_false(self):
        """Test boolean false mapping for probabilistic_sampling."""
        config = {"probabilistic_sampling": False}
        sample = subsample.Sample("test_sample", config, {})

        assert "--no-probabilistic-sampling" in sample.filter_args
        assert sample.filter_args["--no-probabilistic-sampling"] is None
        assert "--probabilistic-sampling" not in sample.filter_args

    def test_list_mapping(self):
        """Test mapping of list values to dict values."""
        config = {
            "group_by": ["region", "country", "year"],
            "query_columns": ["date", "region"],
            "exclude_where": ["region=Antarctica", "country=Unknown"]
        }

        sample = subsample.Sample("test_sample", config, {})

        # Check that list values are stored as lists in the dict
        assert "--group-by" in sample.filter_args
        assert sample.filter_args["--group-by"] == ["region", "country", "year"]

        assert "--query-columns" in sample.filter_args
        assert sample.filter_args["--query-columns"] == ["date", "region"]

        assert "--exclude-where" in sample.filter_args
        assert sample.filter_args["--exclude-where"] == ["region=Antarctica", "country=Unknown"]
