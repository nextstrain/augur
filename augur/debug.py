"""
Debug flags and utilities.

See also :envvar:`AUGUR_DEBUG`.
"""
from os import environ

DEBUGGING = bool(environ.get("AUGUR_DEBUG"))
