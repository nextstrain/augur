"""
Helpers for YAML-based configuration files.
"""
from pathlib import Path
from textwrap import dedent
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
