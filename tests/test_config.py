import argparse

import pytest

from augur.config import RecordConfigStore, RecordConfigStoreConst, apply_config
from augur.errors import AugurError


CONFIG_OPTION_TYPES = {
    "flag": bool,
    "no_flag": (bool, "flag", True),
    "two_words": bool,
    "count": int,
    "names": list,
    "numbers": list[int],
}

def parse_args(config_file, extra_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    parser.add_argument("--flag", action=RecordConfigStoreConst, const=True, default=False)
    parser.add_argument("--no-flag", action=RecordConfigStoreConst, const=False, dest="flag")
    parser.add_argument("--two-words", action=RecordConfigStoreConst, const=True, default=False)
    parser.add_argument("--count", action=RecordConfigStore, type=int, default=1)
    parser.add_argument("--names", action=RecordConfigStore, nargs="+")
    parser.add_argument("--numbers", action=RecordConfigStore, nargs="+", type=int)

    return parser.parse_args([
        "--config", str(config_file),
        *(extra_args or []),
    ])


def test_apply_config(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("""\
flag: false
count: 5
names:
  - one
  - two
numbers:
  - 1
  - 2
""")

    args = apply_config(parse_args(config), CONFIG_OPTION_TYPES, "test-command")

    assert args.flag is False
    assert args.count == 5
    assert args.names == ["one", "two"]
    assert args.numbers == [1, 2]


def test_apply_config_rejects_cli_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("count: 5\n")

    args = parse_args(config, ["--count", "6"])

    with pytest.raises(AugurError, match="--count"):
        apply_config(args, CONFIG_OPTION_TYPES, "test-command")


def test_apply_config_allows_cli_options_missing_from_config(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("count: 5\n")

    args = apply_config(
        parse_args(config, ["--flag"]),
        CONFIG_OPTION_TYPES,
        "test-command",
    )

    assert args.flag is True
    assert args.count == 5


def test_apply_config_rejects_unknown_config_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("unknown: true\n")

    args = parse_args(config)

    with pytest.raises(AugurError, match="Invalid test-command config"):
        apply_config(args, CONFIG_OPTION_TYPES, "test-command")


def test_apply_config_supports_underscored_config_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("two_words: true\n")

    args = apply_config(parse_args(config), CONFIG_OPTION_TYPES, "test-command")

    assert args.two_words is True


def test_apply_config_supports_inverted_config_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("no_flag: true\n")

    args = apply_config(parse_args(config), CONFIG_OPTION_TYPES, "test-command")

    assert args.flag is False


def test_apply_config_rejects_cli_options_for_inverted_config_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("no_flag: true\n")

    args = parse_args(config, ["--flag"])

    with pytest.raises(AugurError, match="--flag"):
        apply_config(args, CONFIG_OPTION_TYPES, "test-command")


def test_apply_config_rejects_inverted_cli_options_for_inverted_config_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("no_flag: true\n")

    args = parse_args(config, ["--no-flag"])

    with pytest.raises(AugurError, match="--no-flag"):
        apply_config(args, CONFIG_OPTION_TYPES, "test-command")


def test_apply_config_rejects_duplicate_config_option_destinations(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("""\
flag: true
no_flag: true
""")

    args = parse_args(config)

    with pytest.raises(AugurError, match="multiple config keys map to the same option"):
        apply_config(args, CONFIG_OPTION_TYPES, "test-command")


def test_apply_config_rejects_hyphenated_config_options(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("two-words: true\n")

    args = parse_args(config)

    with pytest.raises(AugurError, match="Invalid test-command config"):
        apply_config(args, CONFIG_OPTION_TYPES, "test-command")
