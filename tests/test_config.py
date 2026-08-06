import pytest
import yaml
from pathlib import Path
from textwrap import dedent
from augur import make_parser


def write_config_file(tmp_path: Path, config: dict) -> Path:
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        yaml.safe_dump(config, f)
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
    config_file = write_config_file(tmp_path, {"invalid": 0})

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

          invalid
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
