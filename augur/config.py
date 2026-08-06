"""
Helpers for YAML-based configuration files.
"""
from argparse import ArgumentParser
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict, List, TypedDict
from augur.argparse_ import ExtendOverwriteDefault, config_key_to_cli_option
from augur.errors import AugurError
from augur.io.print import indented_list

def resolve_filepath(
    path: Path,
    search_paths: list[Path],
) -> Path:
    """
    Resolve a filepath by searching through multiple directories.

    Parameters
    ----------
    path
        The filepath to resolve. May be either an absolute path or a path
        relative to one of the directories in ``search_paths``.
    search_paths
        Directories to search, in order, when ``path`` is relative. Ignored when
        ``path`` is absolute.

    Examples
    --------

    If the path is already absolute, verify it exists and return it.

    >>> import tempfile
    >>> tmpdir1 = Path(tempfile.mkdtemp()).resolve()
    >>> tmpdir2 = Path(tempfile.mkdtemp()).resolve()
    >>> absolute_path = tmpdir1 / "file.txt"
    >>> with open(absolute_path, "w") as f: _ = f.write("test")
    >>> resolve_filepath(absolute_path, []) == absolute_path
    True

    Otherwise, try resolving it relative to each directory in search_paths, in order.
    Return the first path that exists.

    >>> with open(tmpdir2 / "file.txt", "w") as f: _ = f.write("test")
    >>> result = resolve_filepath(Path("file.txt"), [tmpdir1, tmpdir2])
    >>> result == tmpdir1 / "file.txt"
    True

    If an absolute path doesn't exist, raise an error.

    >>> resolve_filepath(Path("/nonexistent/file.txt"), [tmpdir1, tmpdir2])
    Traceback (most recent call last):
      ...
    augur.errors.AugurError: File '/nonexistent/file.txt' does not exist.

    If the relative path doesn't exist anywhere, raise an error.

    >>> resolve_filepath(Path("nonexistent.txt"), [tmpdir1, tmpdir2]) # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    augur.errors.AugurError: File 'nonexistent.txt' not resolvable from any of the following paths:
    <BLANKLINE>
      ...
    """
    # Absolute path
    if path.is_absolute():
        if not path.exists():
            raise AugurError(f"File {str(path)!r} does not exist.")
        return path

    # Relative path
    for search_path in search_paths:
        resolved_path = (search_path / path).resolve()
        if resolved_path.exists():
            return resolved_path

    # File not found
    raise AugurError(dedent(f"""\
        File {str(path)!r} not resolvable from any of the following paths:

          {indented_list([str(p) for p in search_paths], '        ' + '  ')}"""))


class OptionInfo(TypedDict, total=False):
    """
    Structure definition for information about an option.
    """
    type: str
    description: str
    is_file: bool
    default: Any
    choices: List[Any]
    required: bool
    mutually_exclusive_group: str
    negates: str

    # only for CLI
    additional_cli_flags: List[str]

    # only for config schema
    oneOf: List[Dict[str, Any]]


CommandOptions = Dict[str, OptionInfo]
"""Dictionary of option name to information."""


CLI_TYPE_MAP: Dict[str, type] = {
    "string": str,
    "integer": int,
    "number": float,
    "boolean": bool,
    "list of strings": str,
    "list of integers": int,
}


SCHEMA_TYPE_MAP: Dict[str, Dict[str, Any]] = {
    "string": {"type": "string"},
    "integer": {"type": "integer"},
    "number": {"type": "number"},
    "boolean": {"type": "boolean"},
    "list of strings": {
        "oneOf": [
            {"type": "string"},
            {
                "type": "array",
                "items": {"type": "string"}
            }
        ]
    },
    "list of integers": {
        "type": "array",
        "items": {"type": "integer"}
    },
}


def add_arguments(parser: ArgumentParser, options: CommandOptions):
    """
    Register CLI options defined in an options dictionary onto the given
    argparse parser.
    """
    # Handle mutually exclusive groups using a mapping of name to group object.
    groups = {}

    for key, spec in options.items():
        # Set mutually exclusive group
        group = None
        if "mutually_exclusive_group" in spec:
            group_name = spec["mutually_exclusive_group"]
            if group_name not in groups:
                groups[group_name] = parser.add_mutually_exclusive_group()
            group = groups[group_name]

        args = [config_key_to_cli_option(key), *spec.get("additional_cli_flags", [])]
        kwargs: Dict[str, Any] = {}

        # Handle file options
        if spec.get("is_file"):
            kwargs["type"] = str
            kwargs["metavar"] = "FILE"

        # Handle special type mappings and booleans
        if spec.get("type") in ("list of strings", "list of integers"):
            kwargs["type"] = CLI_TYPE_MAP[spec["type"]]
            kwargs["nargs"] = "+"
            kwargs["action"] = ExtendOverwriteDefault
        elif spec.get("type") == "boolean":
            if "negates" in spec:
                kwargs["action"] = "store_false"
                kwargs["dest"] = spec["negates"]
            else:
                kwargs["action"] = "store_true"
        elif "type" in spec:
            kwargs["type"] = CLI_TYPE_MAP[spec["type"]]

        # Handle other fields
        if "description" in spec:
            kwargs["help"] = spec["description"]
        if "choices" in spec:
            kwargs["choices"] = spec["choices"]
        if "default" in spec:
            kwargs["default"] = spec["default"]
        if "required" in spec:
            kwargs["required"] = spec["required"]

        (group or parser).add_argument(*args, **kwargs)
