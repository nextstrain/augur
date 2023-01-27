import argparse
import datetime
from textwrap import dedent
import isodate
import pandas as pd
import re
import treetime.utils
from .errors import AugurError

from augur.util_support.date_disambiguator import DateDisambiguator

SUPPORTED_DATE_HELP_TEXT = dedent("""\
    1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
    2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
    3. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
""")

def numeric_date(date):
    """
    Converts the given *date* to a :py:class:`float`.

    Parameters
    ----------
    date
        Date in any of the supported formats.


    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04")
    2020.42486...
    >>> import datetime, isodate, treetime
    >>> numeric_date("1W") == treetime.utils.numeric_date(datetime.date.today() - isodate.parse_duration("P1W"))
    True
    """
    # Absolute date as a datetime.date.
    if isinstance(date, datetime.date):
        # Use a treetime utility function to convert the datetime.date to a
        # numeric representation.
        return treetime.utils.numeric_date(date)

    # Absolute date in numeric format.
    try:
        return float(date)
    except ValueError:
        pass

    # Absolute date in ISO 8601 date format.
    try:
        return treetime.utils.numeric_date(datetime.date(*map(int, date.split("-", 2))))
    except ValueError:
        pass

    # Backwards-looking relative date in ISO 8601 duration format.
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
    """Get the numeric date from any supported date format.

    This function is intended to be used as the `type` parameter in `argparse.ArgumentParser.add_argument()`

    This raises an ArgumentTypeError, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    """
    try:
        return numeric_date(date)
    except ValueError as e:
        raise argparse.ArgumentTypeError(str(e)) from e

def ambiguous_date_to_date_range(uncertain_date, fmt, min_max_year=None):
    return DateDisambiguator(uncertain_date, fmt=fmt, min_max_year=min_max_year).range()

def is_date_ambiguous(date, ambiguous_by="any"):
    """
    Returns whether a given date string in the format of YYYY-MM-DD is ambiguous by a given part of the date (e.g., day, month, year, or any parts).

    Parameters
    ----------
    date : str
        Date string in the format of YYYY-MM-DD
    ambiguous_by : str
        Field of the date string to test for ambiguity ("day", "month", "year", "any")
    """
    date_components = date.split('-', 2)

    if len(date_components) == 3:
        year, month, day = date_components
    elif len(date_components) == 2:
        year, month = date_components
        day = "XX"
    else:
        year = date_components[0]
        month = "XX"
        day = "XX"

    # Determine ambiguity hierarchically such that, for example, an ambiguous
    # month implicates an ambiguous day even when day information is available.
    return any((
        "X" in year,
        "X" in month and ambiguous_by in ("any", "month", "day"),
        "X" in day and ambiguous_by in ("any", "day")
    ))

def get_numerical_date_from_value(value, fmt=None, min_max_year=None):
    value = str(value)
    if re.match(r'^-*\d+\.\d+$', value):
        # numeric date which can be negative
        return float(value)
    if value.isnumeric():
        # year-only date is ambiguous
        value = fmt.replace('%Y', value).replace('%m', 'XX').replace('%d', 'XX')
    if 'XX' in value:
        ambig_date = ambiguous_date_to_date_range(value, fmt, min_max_year)
        if ambig_date is None or None in ambig_date:
            return [None, None] #don't send to numeric_date or will be set to today
        return [treetime.utils.numeric_date(d) for d in ambig_date]
    try:
        return treetime.utils.numeric_date(datetime.datetime.strptime(value, fmt))
    except:
        return None

def get_numerical_dates(metadata:pd.DataFrame, name_col = None, date_col='date', fmt=None, min_max_year=None):
    if not isinstance(metadata, pd.DataFrame):
        raise AugurError("Metadata should be a pandas.DataFrame.")
    if fmt:
        strains = metadata.index.values
        dates = metadata[date_col].apply(
            lambda date: get_numerical_date_from_value(
                date,
                fmt,
                min_max_year
            )
        ).values
    else:
        strains = metadata.index.values
        dates = metadata[date_col].astype(float)
    return dict(zip(strains, dates))

def get_iso_year_week(year, month, day):
    return datetime.date(year, month, day).isocalendar()[:2]
