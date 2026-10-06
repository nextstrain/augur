"""
Helpers for YAML-based configuration files.
"""
import os
from pathlib import Path
from ruamel.yaml import YAML, YAMLError, constructor
from textwrap import dedent
from typing import Any, Optional
from augur.errors import AugurError
from augur.io.print import indented_list
from augur.validate import load_augur_json_schema, ValidateError

COMMAND_SCHEMAS = {
    "align": "v1",
    "ancestral": "v1",
    "frequencies": "v1",
    "mask": "v1",
    "refine": "v1",
    "traits": "v1",
    "translate": "v1",
    "tree": "v1",
}


def get_referenced_files(config_file: str | Path) -> set[str]:
    """Get the files referenced in an augur config file.

    Extracts and resolves all filepath values referenced in the config.

    Parameters
    ----------
    config_file
        Path to the config file.

    Returns
    -------
    set
        Resolved filepaths
    """
    config = parse_config(config_file)

    schema_ref = config.get("$schema")
    if not schema_ref:
        raise AugurError(f"The configuration file {str(config_file)!r} does not specify a '$schema'.")

    try:
        schema_validator = load_augur_json_schema(str(schema_ref))
    except (FileNotFoundError, ValidateError) as e:
        raise AugurError(f"Schema {schema_ref!r} not found: {e}") from e

    # Resolve filepaths.
    search_paths = get_search_paths(config_file)
    config, filepaths = resolve_filepaths(config, search_paths, schema_validator.schema)

    return set(filepaths)


def parse_config(filename: str | Path) -> dict[str, Any]:
    # Create a custom YAML constructor to treat timestamps as strings.
    class CustomConstructor(constructor.SafeConstructor):
        pass
    def string_constructor(loader, node):
        return loader.construct_scalar(node)
    CustomConstructor.add_constructor('tag:yaml.org,2002:timestamp', string_constructor)

    yaml = YAML(typ="safe")
    yaml.Constructor = CustomConstructor

    with open(filename) as f:
        try:
            config = yaml.load(f)
        except YAMLError as e:
            raise AugurError(f"The configuration file {filename!r} is not valid YAML.\n" + str(e)) from e

    return config


def get_search_paths(config_file: str | Path) -> list[Path]:
    """
    Returns the paths to search for relative filepaths in config.
    """
    default = [
        Path(config_file).parent,
        Path.cwd(),
    ]

    from_env = os.environ.get('AUGUR_SEARCH_PATHS')

    if from_env:
        return [
            *(Path(p) for p in from_env.split(':')),
            *default,
        ]

    return default


def resolve_filepaths(
    config: dict[str, Any],
    search_paths: list[Path],
    schema: dict[str, Any],
    root_schema: Optional[dict[str, Any]] = None,
) -> tuple[dict[str, Any], list[str]]:
    """
    Resolve filepaths in config.

    Recursively walks the config alongside the schema to determine which fields
    contain filepaths, resolves them, and collects the resolved filepaths.
    """
    if root_schema is None:
        root_schema = schema

    filepaths = []

    # Get properties schema for current section
    properties = schema.get("properties", {})
    pattern_properties = schema.get("patternProperties", {})

    for key, value in config.items():
        if key == "$schema" or key.startswith("_"):
            continue
        prop_schema = properties.get(key)

        if not prop_schema and pattern_properties:
            # Use first pattern property schema (for dynamic keys like samples)
            prop_schema = next(iter(pattern_properties.values()))

        # Get referenced property schema
        if ref := prop_schema.get("$ref"):
            prop_schema = _get_referenced_schema(ref, root_schema)
        elif "oneOf" in prop_schema and isinstance(value, dict):
            _, prop_schema = best_matching_variant(prop_schema["oneOf"], value, root_schema)

        # Resolve filepath
        if _is_filepath(prop_schema):
            if isinstance(value, list):
                config[key] = [str(resolve_filepath(Path(v), search_paths)) for v in value]
                filepaths.extend(config[key])
            elif isinstance(value, str):
                config[key] = str(resolve_filepath(Path(value), search_paths))
                filepaths.append(config[key])

        # Recurse into config section
        elif isinstance(value, dict):
            config[key], downstream_filepaths = resolve_filepaths(
                value, search_paths, prop_schema, root_schema
            )
            filepaths.extend(downstream_filepaths)

    return config, filepaths


def best_matching_variant(
    variants: list[dict[str, Any]],
    value: dict[str, Any],
    root_schema: dict[str, Any],
) -> tuple[str, dict[str, Any]]:
    """
    Given a oneOf list of schema variants (typically $ref entries), resolve each
    and return the one whose properties best match the keys in *value*, as a
    (name, schema) tuple.

    The name is the last component of the matched ``$ref`` (e.g.
    ``"filterSampleProperties"``).
    """
    value_keys = set(value.keys())
    best_schema = None
    best_ref = None
    best_overlap = -1
    for variant in variants:
        ref = variant.get("$ref")
        if ref:
            resolved = _get_referenced_schema(ref, root_schema)
        else:
            resolved = variant
        props = set(resolved.get("properties", {}).keys())
        overlap = len(value_keys & props)
        if overlap > best_overlap:
            best_overlap = overlap
            best_schema = resolved
            best_ref = ref
    if not best_schema:
        raise AugurError("Couldn't match oneOf schema for config dict")
    return (best_ref.rsplit("/", 1)[-1], best_schema)


def _get_referenced_schema(
    ref: str,
    root_schema: dict[str, Any],
) -> dict[str, Any]:
    """
    Resolve a JSON schema reference. Example: '#/$defs/filterSampleProperties'
    """
    keys = ref.lstrip("#/").split("/")
    schema = root_schema
    for key in keys:
        schema = schema[key]
    return schema


def _is_filepath(prop_schema: dict[str, Any]) -> bool:
    """
    Check if the property schema declares it is a filepath.
    """
    # Direct 'format: filepath'
    if prop_schema.get("format") == "filepath":
        return True

    # Check oneOf variants for 'format: filepath'
    if "oneOf" in prop_schema:
        for variant in prop_schema["oneOf"]:
            if variant.get("format") == "filepath":
                return True

    return False


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
