import argparse
import datetime
import re
import isodate
import pandas as pd
import re
import treetime.utils
from textwrap import dedent
from functools import lru_cache
from typing import Any


class InvalidDateFormat(ValueError):
    pass


from augur.util_support.date_disambiguator import DateDisambiguator

SUPPORTED_DATE_HELP_TEXT = dedent("""\
    1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
    2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
    3. an ambiguous date in ISO 8601-like format (e.g. '2020-06-XX', '2020-XX-XX') or
    4. an incomplete date in ISO 8601-like format (e.g. '2020-06', '2020') or
    5. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
""")

# float, negative ok
# note: year-only is treated as incomplete ambiguous and must be non-negative (see RE_YEAR_ONLY)
RE_NUMERIC_DATE = re.compile(r'^-*\d+\.\d+$')

# complete ISO 8601 date
# e.g. 2018-03-25
RE_ISO_8601_DATE = re.compile(r'^\d{4}-\d{2}-\d{2}$')

# complete ambiguous ISO 8601 date
# e.g. 2018-03-XX
RE_AMBIGUOUS_ISO_8601_DATE = re.compile(r'^[\dX]{4}-[\dX]{2}-[\dX]{2}$')

# incomplete ambiguous ISO 8601 date (missing day)
# e.g. 2018-03
RE_AMBIGUOUS_ISO_8601_DATE_YEAR_MONTH = re.compile(r'^[\dX]{4}-[\dX]{2}$')

# incomplete ambiguous ISO 8601 date (missing month and day)
# e.g. 2018
# and other positive ints
# e.g. 1, 123, 12345
RE_YEAR_ONLY = re.compile(r'^[1-9][\dX]*$')

# also support relative dates (ISO 8601 durations) - see any_to_numeric()


CACHE_SIZE = 8192
# Some functions below use @lru_cache to minimize redundant operations on
# large datasets that are likely to have multiple entries with the same date value.


def iso_to_numeric(date_in:str, ambiguity_resolver:str=None):
    """Convert ISO 8601 date string to numeric, resolving any ambiguity detected by explicit 'X' characters or missing date parts.
    Parameters
    ----------
    date_in
        Date string in ISO 8601 format.
    ambiguity_resolver
        None: Assume the given ISO 8601 date string is an exact date.
        'min': Resolve to minimum of ambiguous range.
        'max': Resolve to maximum of ambiguous range.
    >>> round(iso_to_numeric('2018-01-01'), 3)
    2018.001
    >>> round(iso_to_numeric('2018-03-25'), 2)
    2018.23
    >>> round(iso_to_numeric('2018-03', ambiguity_resolver='min'), 3)
    2018.163
    >>> round(iso_to_numeric('2018-03', ambiguity_resolver='max'), 3)
    2018.245
    >>> round(iso_to_numeric('2018', ambiguity_resolver='min'), 3)
    2018.001
    >>> round(iso_to_numeric('2018', ambiguity_resolver='max'), 3)
    2018.999
    >>> iso_to_numeric('2018-03-XX')
    Traceback (most recent call last):
      ...
    ValueError: invalid literal for int() with base 10: 'XX'
    >>> round(iso_to_numeric('2018-03-XX', ambiguity_resolver='min'), 3)
    2018.163
    >>> round(iso_to_numeric('2018-03-XX', ambiguity_resolver='max'), 3)
    2018.245
    >>> round(iso_to_numeric('2018-XX-25', ambiguity_resolver='min'), 3)
    2018.067
    >>> round(iso_to_numeric('2018-XX-25', ambiguity_resolver='max'), 3)
    2018.982
    >>> round(iso_to_numeric('201X-XX-XX', ambiguity_resolver='min'), 3)
    2010.001
    >>> round(iso_to_numeric('201X-XX-XX', ambiguity_resolver='max'), 3)
    2019.999
    """
    date_parts = date_in.split('-', maxsplit=2)
    # TODO: resolve partial month/day ambiguity eg. 2018-1X-XX, 2018-10-3X
    if ambiguity_resolver is None:
        year = int(date_parts[0])
        month = int(date_parts[1])
        day = int(date_parts[2])
    elif ambiguity_resolver == 'min':
        year = int(date_parts[0].replace('X', '0'))
        month = int(date_parts[1]) if len(date_parts) > 1 and date_parts[1].isnumeric() else 1
        day = int(date_parts[2]) if len(date_parts) > 2 and date_parts[2].isnumeric() else 1
    elif ambiguity_resolver == 'max':
        year = int(date_parts[0].replace('X', '9'))
        month = int(date_parts[1]) if len(date_parts) > 1 and date_parts[1].isnumeric() else 12
        if len(date_parts) == 3 and date_parts[2].isnumeric():
            day = int(date_parts[2])
        else:
            if month in {1,3,5,7,8,10,12}:
                day = 31
            elif month == 2:
                day = 28
            else:
                day = 30
    try:
        return date_to_numeric_capped(datetime.date(year, month, day))
    except ValueError:
        # catches month/day out of bounds errors
        raise InvalidDateFormat


def relative_iso_to_numeric(backwards_duration_str:str, from_date:datetime.date=None):
    """Compute numeric date from a ISO 8601 duration string relative to a specified date.
    Parameters
    ----------
    backwards_duration_str
        ISO 8601 duration string specifying the duration to go back. The 'P' duration designator is optional.
    from_date
        The date to go back from. Default is current date.
    >>> round(relative_iso_to_numeric('5D', from_date=datetime.date(2018, 3, 25)), 3)
    2018.215
    >>> round(relative_iso_to_numeric('P5D', from_date=datetime.date(2018, 3, 25)), 3)
    2018.215
    >>> round(relative_iso_to_numeric('5W', from_date=datetime.date(2018, 3, 25)), 3)
    2018.133
    >>> round(relative_iso_to_numeric('5Y', from_date=datetime.date(2018, 3, 25)), 3)
    2013.229
    """
    if from_date is None:
        from_date = datetime.date.today()
    if not backwards_duration_str.startswith('P'):
        backwards_duration_str = 'P'+backwards_duration_str
    return treetime.utils.numeric_date(from_date - isodate.parse_duration(backwards_duration_str))


def any_to_numeric(date_in:Any, ambiguity_resolver:str=None):
    """Return numeric date if date is in a supported format.
    For ambiguous ISO 8601 dates, resolve to either minimum or maximum possible value.
    Parameters
    ----------
    date_in
        Date string in any of the supported formats.
    ambiguity_resolver
        None: Assume the given ISO 8601 date string is an exact date.
        'min': Resolve to minimum of ambiguous range.
        'max': Resolve to maximum of ambiguous range.
    >>> round(any_to_numeric(2018, ambiguity_resolver='min'), 3)
    2018.001
    >>> round(any_to_numeric(2018, ambiguity_resolver='max'), 3)
    2018.999
    >>> round(any_to_numeric(2018.0, ambiguity_resolver='max'), 3)
    2018.0
    """
    date_in = str(date_in)

    # Absolute date in numeric format.
    if RE_NUMERIC_DATE.match(date_in):
        return float(date_in)

    # Absolute date in potentially incomplete/ambiguous ISO 8601 date format.
    if (RE_ISO_8601_DATE.match(date_in) or
        RE_AMBIGUOUS_ISO_8601_DATE.match(date_in) or
        RE_AMBIGUOUS_ISO_8601_DATE_YEAR_MONTH.match(date_in) or
        RE_YEAR_ONLY.match(date_in)
        ):
        return iso_to_numeric(date_in, ambiguity_resolver)

    # Relative date in ISO 8601 duration format.
    # No regex for this (it is complex), just try evaluating last and
    # let any expected errors pass to raise the general-purpose InvalidDateFormat.
    try:
        return relative_iso_to_numeric(date_in)
    except (ValueError, isodate.ISO8601Error):
        pass

    raise InvalidDateFormat(f"""Unable to determine date from '{date_in}'. Ensure it is in one of the supported formats:\n{SUPPORTED_DATE_HELP_TEXT}""")


def any_to_numeric_type_min(date_in:Any):
    """Get the numeric date from any supported date format, taking the minimum possible value if ambiguous.
    This function is intended to be used as the `type` parameter in `argparse.ArgumentParser.add_argument()`
    This raises an ArgumentTypeError from InvalidDateFormat exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    >>> round(any_to_numeric_type_min(2018), 3)
    2018.001
    """
    try:
        return any_to_numeric(date_in, ambiguity_resolver='min')
    except InvalidDateFormat as e:
        raise argparse.ArgumentTypeError(str(e)) from e


def any_to_numeric_type_max(date_in:Any):
    """Get the numeric date from any supported date format, taking the maximum possible value if ambiguous.
    This function is intended to be used as the `type` parameter in `argparse.ArgumentParser.add_argument()`
    This raises an ArgumentTypeError from InvalidDateFormat exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    >>> round(any_to_numeric_type_max(2018), 3)
    2018.999
    """
    try:
        return any_to_numeric(date_in, ambiguity_resolver='max')
    except InvalidDateFormat as e:
        raise argparse.ArgumentTypeError(str(e)) from e


### date_to_numeric logic ###
# copied from treetime.utils.numeric_date
# simplified+cached for speed

from calendar import isleap
def date_to_numeric(date:datetime.date):
    """Return the numeric date representation of a datetime.date."""
    days_in_year = 366 if isleap(date.year) else 365
    return date.year + (date.timetuple().tm_yday-0.5) / days_in_year


@lru_cache(maxsize=CACHE_SIZE)
def date_to_numeric_capped(date:datetime.date, max_numeric:float=None):
    """Return the numeric date representation of a datetime.date, capped at a maximum numeric value.
    Parameters
    ----------
    date
        Date to convert.
    max_numeric
        Maximum numeric value to return. Default is current date.
    >>> round(date_to_numeric_capped(datetime.date(2017, 1, 1), max_numeric=2018.0), 3)
    2017.001
    >>> round(date_to_numeric_capped(datetime.date(2019, 1, 1), max_numeric=2018.0), 3)
    2018.0
    """
    if max_numeric is None:
        max_numeric = date_to_numeric(datetime.date.today())
    date_numeric = date_to_numeric(date)
    if date_numeric > max_numeric:
        date_numeric = max_numeric
    return date_numeric


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
        return {}
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
