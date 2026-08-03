"""
Helpers for YAML-based configuration files.
"""
import os
from pathlib import Path
from textwrap import dedent
from typing import Any, Dict, List, Optional, Tuple
from referencing import Resource
from augur.errors import AugurError
from augur.io.print import indented_list, print_err
from augur.validate import create_local_schema_registry


def resolve_filepaths(
    config: Dict[str, Any],
    search_paths: List[Path],
    schema: Dict[str, Any],
    root_schema: Optional[Dict[str, Any]] = None,
) -> Tuple[Dict[str, Any], List[str]]:
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
        if key.startswith("_"):
            continue
        prop_schema = properties.get(key)

        if not prop_schema and pattern_properties:
            # Use first pattern property schema (for dynamic keys like samples)
            prop_schema = next(iter(pattern_properties.values()))

        # Get referenced property schema
        if ref := prop_schema.get("$ref"):
            prop_schema = _get_referenced_schema(ref, root_schema)
            # Update root for external schemas.
            # References within these are relative to itself rather than the outer schema.
            if prop_schema.get("$id"):
                root_schema = prop_schema
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


def get_search_paths(
    config_file: str,
    from_cli: List[str],
) -> List[Path]:
    """
    Returns the paths to search for relative filepaths in config.
    """
    default = [
        Path(config_file).parent,
        Path.cwd(),
    ]

    from_env = os.environ.get('AUGUR_SEARCH_PATHS')

    if from_cli:
        if from_env:
            print_err(dedent(f"""\
                WARNING: Both the command line argument --search-paths
                and the environment variable AUGUR_SEARCH_PATHS are set.
                Only the command line argument will be used."""))
        return [
            *(Path(p) for p in from_cli),
            *default,
        ]

    if from_env:
        return [
            *(Path(p) for p in from_env.split(':')),
            *default,
        ]

    return default


def best_matching_variant(
    variants: List[Dict[str, Any]],
    value: Dict[str, Any],
    root_schema: Dict[str, Any],
) -> Tuple[str, Dict[str, Any]]:
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
    name = best_ref.rsplit("/", 1)[-1] if best_ref else ""
    return (name, best_schema)


def _get_referenced_schema(
    ref: str,
    root_schema: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Resolve a JSON schema reference.
    Examples:
        '#/$defs/filterSampleProperties'
        'https://nextstrain.org/schemas/augur/subsample-config/v1'
        'https://nextstrain.org/schemas/augur/subsample-config/v1#/$defs/filterSampleProperties'
    """
    root_id = root_schema.get("$id", "")
    registry = create_local_schema_registry()
    registry = registry.with_resource(root_id, Resource.from_contents(root_schema))
    resolver = registry.resolver(base_uri=root_id)
    return resolver.lookup(ref).contents


def _is_filepath(prop_schema: Dict[str, Any]) -> bool:
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
