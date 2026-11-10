"""
Helpers for YAML-based configuration files.
"""
import os
from pathlib import Path
from ruamel.yaml import YAML, YAMLError, constructor
from textwrap import dedent
from typing import Any, Iterator, Optional
from augur.errors import AugurError
from augur.io.print import indented_list
from augur.schema import InvalidSchemaError, load_augur_json_schema


def get_referenced_files(config_file: str | Path, schema=None) -> set[str]:
    """Get the files referenced in an augur config file.

    Extracts and resolves all path values referenced in the config. Directory
    paths are recursively expanded to all files within.

    Parameters
    ----------
    config_file
        Path to the config file.
    schema
        Augur schema URL or file name. Default: infer from $schema

    Returns
    -------
    set
        Resolved filepaths
    """
    config = parse_config(config_file)

    schema_ref = schema or config.get("$schema")
    if not schema_ref:
        raise AugurError(f"The configuration file {str(config_file)!r} does not specify a '$schema'.")

    try:
        schema_validator = load_augur_json_schema(str(schema_ref))
    except FileNotFoundError as e:
        raise AugurError(f"Schema {schema_ref!r} not found: {e}") from e
    except InvalidSchemaError as e:
        raise AugurError(str(e)) from e

    # Resolve paths.
    search_paths = get_search_paths(config_file)
    config, paths = resolve_paths(config, search_paths, schema_validator.schema)

    return {str(path) for path in _expand_directories(paths)}


def _expand_directories(paths: list[Path]) -> Iterator[Path]:
    """
    Replace directory paths with the paths of all files within them.
    """
    for path in paths:
        if path.is_dir():
            # TODO: use Path.walk() once Python 3.12 is our minimum.
            for file in path.rglob("*"):
                if file.is_file():
                    yield file
        else:
            # Most likely normal files, but could also be special files (e.g. /dev/null).
            yield path


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
    Returns the paths to search for relative paths in config.
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


def resolve_paths(
    config: dict[str, Any],
    search_paths: list[Path],
    schema: dict[str, Any],
    root_schema: Optional[dict[str, Any]] = None,
) -> tuple[dict[str, Any], list[Path]]:
    """
    Resolve paths in config.

    Recursively walks the config alongside the schema to determine which fields
    contain paths, resolves them, and collects the resolved paths.
    """
    if root_schema is None:
        root_schema = schema

    paths = []

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

        # Resolve file or directory path.
        if _is_path(prop_schema):
            if isinstance(value, list):
                resolved = [resolve_path(Path(v), search_paths) for v in value]
                config[key] = [str(path) for path in resolved]
                paths.extend(resolved)
            elif isinstance(value, str):
                resolved_path = resolve_path(Path(value), search_paths)
                config[key] = str(resolved_path)
                paths.append(resolved_path)

        # Recurse into config section
        elif isinstance(value, dict):
            config[key], downstream_paths = resolve_paths(
                value, search_paths, prop_schema, root_schema
            )
            paths.extend(downstream_paths)

    return config, paths


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


PATH_FORMATS = {"filepath", "dirpath"}


def _is_path(prop_schema: dict[str, Any]) -> bool:
    """
    Check if the property schema declares it is a file or directory path.
    """
    # Direct 'format: filepath' or 'format: dirpath'
    if prop_schema.get("format") in PATH_FORMATS:
        return True

    # Check oneOf variants for 'format: filepath' or 'format: dirpath'
    if "oneOf" in prop_schema:
        for variant in prop_schema["oneOf"]:
            if variant.get("format") in PATH_FORMATS:
                return True

    return False


def resolve_path(
    path: Path,
    search_paths: list[Path],
) -> Path:
    """
    Resolve a file or directory path by searching through multiple directories.

    Parameters
    ----------
    path
        The path to resolve. May be either an absolute path or a path
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
    >>> resolve_path(absolute_path, []) == absolute_path
    True

    Otherwise, try resolving it relative to each directory in search_paths, in order.
    Return the first path that exists.

    >>> with open(tmpdir2 / "file.txt", "w") as f: _ = f.write("test")
    >>> result = resolve_path(Path("file.txt"), [tmpdir1, tmpdir2])
    >>> result == tmpdir1 / "file.txt"
    True

    Directories are resolved the same way.

    >>> (tmpdir2 / "dir").mkdir()
    >>> resolve_path(Path("dir"), [tmpdir1, tmpdir2]) == tmpdir2 / "dir"
    True

    If an absolute path doesn't exist, raise an error.

    >>> resolve_path(Path("/nonexistent/file.txt"), [tmpdir1, tmpdir2])
    Traceback (most recent call last):
      ...
    augur.errors.AugurError: Path '/nonexistent/file.txt' does not exist.

    If the relative path doesn't exist anywhere, raise an error.

    >>> resolve_path(Path("nonexistent.txt"), [tmpdir1, tmpdir2]) # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    augur.errors.AugurError: Path 'nonexistent.txt' not resolvable from any of the following paths:
    <BLANKLINE>
      ...
    """
    # Absolute path
    if path.is_absolute():
        if not path.exists():
            raise AugurError(f"Path {str(path)!r} does not exist.")
        return path

    # Relative path
    for search_path in search_paths:
        resolved_path = (search_path / path).resolve()
        if resolved_path.exists():
            return resolved_path

    # Path not found
    raise AugurError(dedent(f"""\
        Path {str(path)!r} not resolvable from any of the following paths:

          {indented_list([str(p) for p in search_paths], '        ' + '  ')}"""))


# TODO: remove this re-export after downstream workflows have been updated to
# use resolve_path.
def resolve_filepath(
    path: Path,
    search_paths: list[Path],
) -> Path:
    return resolve_path(path, search_paths)
