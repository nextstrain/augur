import argparse
import datetime
from textwrap import dedent
import isodate
import pandas as pd
import re
from functools import cache
from treetime.utils import numeric_date as tt_numeric_date, datetime_from_numeric
from typing import Any, Dict, Optional, Tuple, Union
from augur.errors import AugurError
from .errors import InvalidDate

from .ambiguous_date import AmbiguousDate

SUPPORTED_DATE_HELP_TEXT = dedent("""\
    1. an Augur-style numeric date with the year as the integer part (e.g.
       2020.42) or
    2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
    3. a backwards-looking relative date in ISO 8601 duration format with
       optional P prefix (e.g. '1W', 'P1W')""")

def date_to_numeric(date: datetime.date) -> float:
    """Wrapper around treetime.utils.numeric_date that ensures a float is returned."""
    value = tt_numeric_date(date)
    if not isinstance(value, float):
        raise Exception("treetime.utils.numeric_date unexpectedly returned a non-float.")
    return value

def numeric_date(date):
    """
    Converts the given *date* string to a :py:class:`float`.

    *date* may be given as:
    1. A string or float (number) with year as the integer part
    2. A string in the YYYY-MM-DD (ISO 8601) syntax
    3. A string representing a relative date (duration before datetime.date.today())

    Examples
    --------
    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04")
    2020.42486...
    >>> import datetime, isodate
    >>> numeric_date("1W") == date_to_numeric(datetime.date.today() - isodate.parse_duration("P1W"))
    True
    """
    # date is a datetime.date
    if isinstance(date, datetime.date):
        return date_to_numeric(date)

    # date is numeric
    try:
        return float(date)
    except ValueError:
        pass

    # date is in YYYY-MM-DD form
    try:
        return date_to_numeric(datetime.date(*map(int, date.split("-", 2))))
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
        return date_to_numeric(datetime.date.today() - isodate.parse_duration(duration_str))
    except (ValueError, isodate.ISO8601Error):
        pass

    raise InvalidDate(date, f"""Ensure it is in one of the supported formats:\n{SUPPORTED_DATE_HELP_TEXT}""")

def numeric_date_type(date):
    """Wraps numeric_date() for argparse usage.

    This raises an ArgumentTypeError from InvalidDateFormat exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    """
    try:
        return numeric_date(date)
    except InvalidDate as error:
        raise argparse.ArgumentTypeError(str(error)) from error

def is_date_ambiguous(date, ambiguous_by):
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
        year = date_components[0] if date_components[0] else 'X'
        month = "XX"
        day = "XX"

    # Determine ambiguity hierarchically such that, for example, an ambiguous
    # month implicates an ambiguous day even when day information is available.
    return any((
        "X" in year,
        "X" in month and ambiguous_by in ("any", "month", "day"),
        "X" in day and ambiguous_by in ("any", "day")
    ))

RE_NUMERIC_DATE = re.compile(r'^-*\d+\.\d+$')
"""
Matches floats (e.g. 2018.0, -2018.0).
Note that a year-only value is treated as incomplete ambiguous and must be
non-negative (see :const:`RE_YEAR_ONLY`).
"""

RE_YEAR_ONLY = re.compile(r'^\d+$')
"""
Matches:

1. Incomplete ambiguous ISO 8601 dates that are missing both the month and day
   parts (e.g. 2018)
2. Other positive integers (e.g. 1, 123, 12345)
"""

RE_YEAR_MONTH_ONLY = re.compile(r'^\d+-\d+$')
"""
Matches:

1. Reduced precision ISO 8601 dates that are missing the day part (e.g. 2018-03)
2. Other dates in <year>-<month> format (e.g. 2018-3)

Note: This matches out of bounds dates such as 2018-13.
Those should be further validated by date conversion functions.
"""

RE_YEAR_MONTH_DAY = re.compile(r'^\d+-\d+-\d+$')
"""
Matches:

1. ISO 8601 dates (e.g. 2018-03-09)
2. Other dates in <year>-<month>-<day> format (e.g. 2018-3-9)

Note: This matches out of bounds dates such as 2018-03-32.
Those should be further validated by date conversion functions.
"""

RE_AUGUR_AMBIGUOUS_DATE = re.compile(r'.*XX.*')
"""
Matches an Augur-style ambiguous date with 'XX' used to mask unknown parts of the date.
Note that this can support any date format, not just YYYY-MM-DD.
"""

@cache
def get_numerical_date_from_value(value, fmt, min_max_year=None) -> Union[float, Tuple[float, float], None]:
    value = str(value)

    # 1. Check if value is an exact date in the specified format (fmt).

    try:
        return date_to_numeric(datetime.datetime.strptime(value, fmt))
    except:
        pass
    
    # 2. Check if value is an ambiguous date in the specified format (fmt).

    if RE_AUGUR_AMBIGUOUS_DATE.match(value):
        start, end = AmbiguousDate(value, fmt=fmt).range(min_max_year=min_max_year)
        return (date_to_numeric(start), date_to_numeric(end))

    # 3. Check formats that are always supported.

    if RE_NUMERIC_DATE.match(value):
        return float(value)

    if RE_YEAR_ONLY.match(value):
        start, end = AmbiguousDate(f"{value}-XX-XX", fmt="%Y-%m-%d").range(min_max_year=min_max_year)
        return (date_to_numeric(start), date_to_numeric(end))

    if RE_YEAR_MONTH_ONLY.match(value):
        start, end = AmbiguousDate(f"{value}-XX", fmt="%Y-%m-%d").range(min_max_year=min_max_year)
        return (date_to_numeric(start), date_to_numeric(end))

    if RE_YEAR_MONTH_DAY.match(value):
        try:
            return date_to_numeric(datetime.date(*map(int, value.split('-'))))
        except ValueError as error:
            # Note: This isn't consistent with the out of bounds handling in
            # AmbiguousDate, which auto-converts out of bounds values to the
            # closest in-bound value.
            raise InvalidDate(value, str(error)) from error

    # 4. Return none (silent error) if the date does not match any of the checked formats.

    return None

def get_numerical_dates(
    metadata: pd.DataFrame,
    fmt,
    name_col = None,
    date_col = 'date',
    min_max_year = None,
) -> Dict[str, Union[float, Tuple[float, float], None]]:
    if not isinstance(metadata, pd.DataFrame):
        raise AugurError("Metadata should be a pandas.DataFrame.")

    strains = metadata.index.values
    dates = metadata[date_col].apply(
        lambda date: get_numerical_date_from_value(
            date,
            fmt,
            min_max_year
        )
    ).values

    return dict(zip(strains, dates))

@cache
def get_year_week(year, month, day):
    year, week = datetime.date(year, month, day).isocalendar()[:2]
    return f"{year}-{str(week).zfill(2)}"

@cache
def get_year_month_day(value: Any) -> Tuple[Optional[int], Optional[int], Optional[int]]:
    """
    Extract year, month, and day from a date value.

    Individual date components can be None if unresolvable.
    """
    date_or_range = get_numerical_date_from_value(value, fmt="%Y-%m-%d")

    if date_or_range is None:
        return (None, None, None)

    # Date range
    elif isinstance(date_or_range, tuple):
        start, end = (datetime_from_numeric(date_or_range[0]),
                      datetime_from_numeric(date_or_range[1]))

        # Only use unambiguous date components
        year = start.year if start.year == end.year else None
        month = start.month if start.month == end.month else None
        day = start.day if start.day == end.day else None

        return (year, month, day)

    # Exact date
    else:
        dt = datetime_from_numeric(date_or_range)
        return (dt.year, dt.month, dt.day)
