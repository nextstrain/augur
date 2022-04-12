import argparse
import datetime
from textwrap import dedent
import isodate
import treetime.utils

SUPPORTED_DATE_HELP_TEXT = dedent("""\
    1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
    2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
    3. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
""")

def numeric_date(date):
    """
    Converts the given *date* string to a :py:class:`float`.

    *date* may be given as:
    1. A string or float (number) with year as the integer part
    2. A string in the YYYY-MM-DD (ISO 8601) syntax
    3. A string representing a relative date (duration before datetime.date.today())

    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04")
    2020.42486...
    >>> import datetime, isodate, treetime
    >>> numeric_date("1W") == treetime.utils.numeric_date(datetime.date.today() - isodate.parse_duration("P1W"))
    True
    """
    # date is numeric
    try:
        return float(date)
    except ValueError:
        pass

    # date is in YYYY-MM-DD form
    try:
        return treetime.utils.numeric_date(datetime.date(*map(int, date.split("-", 2))))
    except ValueError:
        pass

    # date is a duration treated as a backwards-looking relative date
    try:
        # make a copy of date for this block
        duration_str = str(date)
        if duration_str.startswith('P'):
            duration_str = duration_str
        else:
            duration_str = 'P'+duration_str
        return treetime.utils.numeric_date(datetime.date.today() - isodate.parse_duration(duration_str))
    except (ValueError, isodate.ISO8601Error):
        pass

    raise ValueError(f"""Unable to determine date from '{date}'. Ensure it is in one of the supported formats:\n{SUPPORTED_DATE_HELP_TEXT}""")

def numeric_date_type(date):
    """Wraps numeric_date() for argparse usage.

    This raises an ArgumentTypeError, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    """
    try:
        return numeric_date(date)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e)) from e
