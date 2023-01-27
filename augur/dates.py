import argparse
import calendar
import datetime
from textwrap import dedent
import isodate
import pandas as pd
import re
import treetime.utils
from .errors import AugurError

from augur.util_support.date_disambiguator import DateDisambiguator


class InvalidDateFormat(ValueError):
    pass


SUPPORTED_DATE_HELP_TEXT = dedent("""\
    1. an Augur-style numeric date with the year as the integer part (e.g. 2020.42) or
    2. a date in ISO 8601 date format (i.e. YYYY-MM-DD) (e.g. '2020-06-04') or
    3. an ambiguous date in ISO 8601-like format (e.g. '2020-06-XX', '2020-XX-XX') or
    4. an incomplete date in ISO 8601-like format (e.g. '2020-06', '2020') or
    5. a backwards-looking relative date in ISO 8601 duration format with optional P prefix (e.g. '1W', 'P1W')
""")

# Matches floats (e.g. 2018.0, -2018.0).
# Note that a year-only value is treated as incomplete ambiguous and must be
# non-negative (see RE_YEAR_ONLY).
RE_NUMERIC_DATE = re.compile(r'^-?[0-9]*\.[0-9]*$')

# Matches complete ISO 8601 dates (e.g. 2018-03-25).
RE_ISO_8601_DATE = re.compile(r'^[0-9]{4}-[0-9]{2}-[0-9]{2}$')

# Matches complete ambiguous ISO 8601 dates (e.g. 2018-03-XX).
RE_AMBIGUOUS_ISO_8601_DATE = re.compile(r'^[0-9X]{4}-[0-9X]{2}-[0-9X]{2}$')

# Matches incomplete ambiguous ISO 8601 dates that are missing the day part (e.g. 2018-03).
RE_YEAR_MONTH_ONLY = re.compile(r'^[0-9X]{4}-[0-9X]{2}$')

# Matches
# 1. Incomplete ambiguous ISO 8601 dates that are missing both the month and day
#    parts (e.g. 2018)
# 2. Other positive integers (e.g. 1, 123, 12345)
RE_YEAR_ONLY = re.compile(r'^0*[1-9][0-9X]*$')

# Relative dates (ISO 8601 durations) are also supported - see numeric_date().

# 'X' followed by a specific digit does not make sense.
RE_INVALID_AMBIGUITY = re.compile(r'.*X[0-9]+.*')


def numeric_date(date, ambiguity_resolver: str = None):
    """
    Converts the given *date* to a :py:class:`float`.

    Parameters
    ----------
    date
        Date in any of the supported formats.

    ambiguity_resolver
        None: Assume the given ISO 8601 date string is an exact date.
        'min': Resolve to minimum of ambiguous range.
        'max': Resolve to maximum of ambiguous range.


    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04")
    2020.42486...
    >>> round(numeric_date(2018, ambiguity_resolver='min'), 3)
    2018.001
    >>> round(numeric_date(2018, ambiguity_resolver='max'), 3)
    2018.999
    >>> round(numeric_date(2018.0, ambiguity_resolver='max'), 3)
    2018.0
    >>> import datetime, isodate, treetime
    >>> numeric_date("1W") == treetime.utils.numeric_date(datetime.date.today() - isodate.parse_duration("P1W"))
    True
    """
    # Absolute date as a datetime.date.
    if isinstance(date, datetime.date):
        # Use a treetime utility function to convert the datetime.date to a
        # numeric representation.
        unbounded_date = treetime.utils.numeric_date(date)
        # Make one more call to the base case.
        return numeric_date(unbounded_date)

    # All other formats are treated as strings.
    date = str(date)

    # Base case: Absolute date in numeric format.
    if RE_NUMERIC_DATE.match(date):
        return float(date)

    # Absolute date in potentially incomplete/ambiguous ISO 8601 date format.
    if (RE_ISO_8601_DATE.match(date) or
        RE_AMBIGUOUS_ISO_8601_DATE.match(date) or
        RE_YEAR_MONTH_ONLY.match(date) or
        RE_YEAR_ONLY.match(date)
        ):
        return iso_to_numeric(date, ambiguity_resolver)

    # Backwards-looking relative date in ISO 8601 duration format.
    # No regex for this - it would be too complex - just try evaluating last and
    # let any errors pass to raise the InvalidDateFormat.
    try:
        return relative_iso_to_numeric(date)
    except (ValueError, isodate.ISO8601Error):
        pass

    raise InvalidDateFormat(f"""Unable to determine date from '{date}'. Ensure it is in one of the supported formats:\n{SUPPORTED_DATE_HELP_TEXT}""")


def iso_to_numeric(date: str, ambiguity_resolver: str = None):
    """Convert ISO 8601 date string to numeric, resolving any ambiguity detected by explicit 'X' characters or missing date parts.

    Parameters
    ----------
    date
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
    >>> iso_to_numeric('2018')
    Traceback (most recent call last):
      ...
    augur.dates.InvalidDateFormat: Ambiguous date provided without specifying how to resolve ambiguity.
    """
    date_parts = date.split('-', maxsplit=2)

    year = date_parts[0]
    month = date_parts[1] if len(date_parts) > 1 else 'XX'
    day = date_parts[2] if len(date_parts) > 2 else 'XX'

    if ambiguity_resolver is None:
        try:
            year = int(year)
            month = int(month)
            day = int(day)
        except ValueError as error:
            raise InvalidDateFormat("Ambiguous date provided without specifying how to resolve ambiguity.")

    else:
        validate_ambiguity(year, month, day)

        if ambiguity_resolver == 'min':
            year = year.replace('X', '0')
            year = int(year)

            month = int(month.replace('X', '0'))
            if month < 1:
                month = 1

            day = int(day.replace('X', '0'))
            if day < 1:
                day = 1

        elif ambiguity_resolver == 'max':
            year = year.replace('X', '9')
            year = int(year)

            max_month = 12
            month = int(month.replace('X', '9'))
            if month > max_month:
                month = max_month

            try:
                max_day = calendar.monthrange(year, month)[1]
            except calendar.IllegalMonthError as error:
                # Month out of bounds is a user error.
                raise InvalidDateFormat(error) from error
            day = int(day.replace('X', '9'))
            if day > max_day:
                day = max_day

    try:
        return numeric_date(datetime.date(year, month, day))
    except ValueError as error:
        # Month/day out of bounds errors are user errors.
        if str(error) == "month must be in 1..12":
            raise InvalidDateFormat(error) from error
        if str(error) == "day is out of range for month":
            raise InvalidDateFormat(error) from error
        raise error


def validate_ambiguity(year: str, month: str, day: str):
    """Validate ambiguity between date parts."""
    # Validate between date parts
    if is_ambiguous(year) and (not is_ambiguous(month) or not is_ambiguous(day)):
        raise InvalidDateFormat("Invalid date: Year contains uncertainty, so month and day must also be uncertain.")
    if is_ambiguous(month) and not is_ambiguous(day):
        raise InvalidDateFormat("Invalid date: Month contains uncertainty, so day must also be uncertain.")

    # Validate within date parts
    for date_part in (year, month, day):
        if RE_INVALID_AMBIGUITY.match(date_part):
            raise InvalidDateFormat("Ambiguity can not be followed by an exact digit.")


def is_ambiguous(date_part: str):
    """Determine if a date part is ambiguous."""
    return 'X' in date_part


def relative_iso_to_numeric(backwards_duration_str: str, from_date: datetime.date = None):
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
    return numeric_date(from_date - isodate.parse_duration(backwards_duration_str))


def numeric_date_type(date):
    """Get the numeric date from any supported date format.

    This function is intended to be used as the `type` parameter in `argparse.ArgumentParser.add_argument()`

    This raises an ArgumentTypeError from InvalidDateFormat exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    """
    try:
        return numeric_date(date)
    except InvalidDateFormat as e:
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
