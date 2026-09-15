import pytest
from pathlib import Path
from ruamel.yaml import YAML
from textwrap import dedent
from augur import make_parser
from augur.config import get_referenced_files
from augur.errors import AugurError


def write_config_file(tmp_path: Path, config: dict) -> Path:
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        YAML(typ="safe").dump(config, f)
    return config_file


def test_config_list(tmp_path):
    """
    Test that list values in --config are correctly mapped.
    """
    config_file = write_config_file(tmp_path, {"year_bounds": [2000, 2020]})

    parser = make_parser()
    args = parser.parse_args([
        "refine",
        "--config", str(config_file),
        "--tree", "tree.nwk",
    ])
    assert args.year_bounds == [2000, 2020]


def test_config_scalar_for_list(tmp_path):
    """
    Test that scalar values in --config can also represent single-item lists.
    """
    config_file = write_config_file(tmp_path, {"root": "ROOT"})

    parser = make_parser()
    args = parser.parse_args([
        "refine",
        "--config", str(config_file),
        "--tree", "tree.nwk",
    ])
    assert args.root == ["ROOT"]


@pytest.mark.parametrize("value", [False, True])
def test_config_boolean(tmp_path, value):
    """
    Test that boolean values in --config are correctly mapped.
    """
    config_file = write_config_file(tmp_path, {"covariance": value})

    parser = make_parser()
    args = parser.parse_args([
        "refine",
        "--config", str(config_file),
        "--tree", "tree.nwk",
    ])
    assert args.covariance is value


def test_config_error_with_cli_same_option(tmp_path, capsys):
    """
    Test that an error is shown when an option used in CLI is also used in
    --config.
    """
    config_file = write_config_file(tmp_path, {"timetree": True})

    parser = make_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([
            "refine",
            "--config", str(config_file),
            "--tree", "tree.nwk",
            "--timetree",
        ])
    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert captured.err == dedent("""\
        ERROR: Options can be specified in either --config or on the CLI, but not both.

        The following option was specified in both:

          timetree (config YAML), --timetree (CLI)
        """)


def test_config_error_with_cli_same_dest(tmp_path, capsys):
    """
    Test that an error is shown when an option used in CLI targets the same
    argparse destination as another option used in --config.
    """
    config_file = write_config_file(tmp_path, {"covariance": True})

    parser = make_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([
            "refine",
            "--config", str(config_file),
            "--tree", "tree.nwk",
            "--no-covariance",
        ])
    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert captured.err == dedent("""\
        ERROR: Options can be specified in either --config or on the CLI, but not both.

        The following option was specified in both:

          covariance (config YAML), --no-covariance (CLI)
        """)


def test_config_error_with_invalid(tmp_path, capsys):
    """
    Test that an error is shown when an invalid option is used in --config.
    """
    config_file = write_config_file(tmp_path, {"no_covariance": True})

    parser = make_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([
            "refine",
            "--tree", "tree.nwk",
            "--config", str(config_file),
        ])
    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert captured.err == dedent("""\
        ERROR: The following invalid option was specified in --config:

          no_covariance
        """)


def test_config_error_with_dashes(tmp_path, capsys):
    """
    Test that an error is shown when a dashed option name is used in --config.
    """
    config_file = write_config_file(tmp_path, {"date-confidence": True})

    parser = make_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([
            "refine",
            "--tree", "tree.nwk",
            "--config", str(config_file),
        ])
    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert captured.err == dedent("""\
        ERROR: The following invalid option was specified in --config:

          date-confidence
        """)


def test_config_error_with_duplicate(tmp_path, capsys):
    """
    Test that an error is shown when a duplicate option is used in --config.
    """
    config_file = tmp_path / "config.yaml"
    config_file.write_text(dedent("""\
        timetree: true
        timetree: false
    """))

    parser = make_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([
            "refine",
            "--tree", "tree.nwk",
            "--config", str(config_file),
        ])
    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert 'found duplicate key "timetree"' in captured.err


def test_get_referenced_files_tree(tmp_path):
    """
    Test get_referenced_files with augur tree config.
    """
    alignment = tmp_path / "alignment.fasta"
    alignment.touch()
    config_file = write_config_file(tmp_path, {
        "$schema": "https://nextstrain.org/schemas/augur/tree-config/v1",
        "alignment": "alignment.fasta",
    })

    assert get_referenced_files(config_file) == {str(alignment.resolve())}


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
        "$schema": "https://nextstrain.org/schemas/augur/tree-config/v1",
        "alignment": "nonexistent.fasta",
    })

    with pytest.raises(AugurError, match="not resolvable from any of the following paths"):
        get_referenced_files(config_file)
