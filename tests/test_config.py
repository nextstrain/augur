import pytest
from pathlib import Path
from ruamel.yaml import YAML
from augur.config import get_referenced_files
from augur.errors import AugurError


def write_config_file(tmp_path: Path, config: dict) -> Path:
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        YAML(typ="safe").dump(config, f)
    return config_file


def test_get_referenced_files_subsample(tmp_path):
    """
    Test get_referenced_files with nested structure in subsample config.
    """
    include_file = tmp_path / "include.txt"
    include_file.touch()
    exclude_file = tmp_path / "exclude.txt"
    exclude_file.touch()

    config_file = write_config_file(tmp_path, {
        "$schema": "https://nextstrain.org/schemas/augur/subsample-config/v1",
        "defaults": {
            "exclude": "exclude.txt",
        },
        "samples": {
            "sample_a": {
                "include": ["include.txt"],
            },
        },
    })

    assert get_referenced_files(config_file) == {str(include_file.resolve()), str(exclude_file.resolve())}


def test_get_referenced_files_missing_schema(tmp_path):
    """
    Test that an error is raised when $schema is missing.
    """
    config_file = write_config_file(tmp_path, {"alignment": "alignment.fasta"})

    with pytest.raises(AugurError, match=r"does not specify a '\$schema'"):
        get_referenced_files(config_file)


def test_get_referenced_files_nonexistent_command(tmp_path):
    """
    Test that an error is raised when the inferred command schema doesn't exist.
    """
    config_file = write_config_file(tmp_path, {
        "$schema": "https://nextstrain.org/schemas/augur/nonexistent-config/v1",
        "alignment": "alignment.fasta",
    })

    with pytest.raises(AugurError, match="not found"):
        get_referenced_files(config_file)


def test_get_referenced_files_missing_file(tmp_path):
    """
    Test that an error is raised when a referenced file cannot be resolved.
    """
    config_file = write_config_file(tmp_path, {
        "$schema": "https://nextstrain.org/schemas/augur/subsample-config/v1",
        "defaults": {
            "exclude": "nonexistent.txt",
        }
    })

    with pytest.raises(AugurError, match="not resolvable from any of the following paths"):
        get_referenced_files(config_file)
