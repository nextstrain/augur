# Definitions in this file are related to date parsing of argparse parameters.

import argparse
from textwrap import dedent
from . import numeric_date
from .errors import InvalidDate


SUPPORTED_DATE_HELP_TEXT = dedent("""\
    1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
    2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
    3. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
""")


def numeric_date_type_min(date):
    """Get the numeric date from any supported date format, taking the minimum possible value if ambiguous.

    This function is intended to be used as the `type` parameter in `argparse.ArgumentParser.add_argument()`

    This raises an ArgumentTypeError from InvalidDate exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513

    >>> round(numeric_date_type_min("2018"), 3)
    2018.001
    """
    try:
        # TODO: support custom formats
        converted_date = numeric_date(date, fmt="%Y-%m-%d", ambiguity_resolver='min')
        if converted_date is None:
            raise InvalidDate(date, f"Ensure it is in one of the supported formats:\n{SUPPORTED_DATE_HELP_TEXT}")
    except InvalidDate as e:
        raise argparse.ArgumentTypeError(str(e)) from e
    return converted_date


def numeric_date_type_max(date):
    """Get the numeric date from any supported date format, taking the maximum possible value if ambiguous.

    This function is intended to be used as the `type` parameter in `argparse.ArgumentParser.add_argument()`

    This raises an ArgumentTypeError from InvalidDate exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513

    >>> round(numeric_date_type_max("2018"), 3)
    2018.999
    """
    try:
        # TODO: support custom formats
        converted_date = numeric_date(date, fmt="%Y-%m-%d", ambiguity_resolver='max')
        if converted_date is None:
            raise InvalidDate(date, f"Ensure it is in one of the supported formats:\n{SUPPORTED_DATE_HELP_TEXT}")
    except InvalidDate as e:
        raise argparse.ArgumentTypeError(str(e)) from e
    return converted_date
