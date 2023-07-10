"""
Resource files.
"""
# Originally copied from nextstrain/cli/resources/__init__.py in the Nextstrain
# CLI project¹.
#
# ¹ <https://github.com/nextstrain/cli/blob/a5dda9c0579ece7acbd8e2c32a4bbe95df7c0bce/nextstrain/cli/resources/__init__.py>

# We gate usage of the stdlib implementation on 3.11 because that's the first
# version with a full adapter for making the new files() / Traversable API
# backwards compatible with importers only providing the original path() /
# ResourceReader API.  The PyPI backport, on the other hand, contains the full
# adapter since 5.3.0, which we declare as our minimum version in setup.py, so
# we use that even on 3.9 and 3.10.
#
# We're using the new API at all because the original one is being deprecated
# and we want to avoid warnings both from the stdlib implementation on 3.11 and
# from the PyPI backport implementation on older Python versions.
#   -trs, 13 Sept 2022
import sys

if sys.version_info >= (3, 11):
    from importlib.resources import files as _files, as_file as _as_file
else:
    from importlib_resources import files as _files, as_file as _as_file

from pathlib import Path
from typing import ContextManager


def as_file(path: str) -> ContextManager[Path]:
    return _as_file(_files(__name__) / path)
