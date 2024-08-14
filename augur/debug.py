"""
Debug flags and utilities.

.. envvar:: AUGUR_DEBUG

    Set to a truthy value (e.g. 1) to print more information about (handled)
    errors.  For example, when this is not set or falsey, stack traces and
    parent exceptions in an exception chain are omitted from handled errors.
"""
from os import environ

DEBUGGING = bool(environ.get("AUGUR_DEBUG"))
