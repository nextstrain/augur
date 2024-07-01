import datetime
import isodate
import pandas as pd
import re
import treetime.utils
from augur.errors import AugurError
from .errors import InvalidDate, InvalidDateMessage

from .ambiguous_date import AmbiguousDate


# Matches floats (e.g. 2018.0, -2018.0).
# Note that a year-only value is treated as incomplete ambiguous and must be
# non-negative (see RE_YEAR_ONLY).
RE_NUMERIC_DATE = re.compile(r'^-?[0-9]*\.[0-9]*$')

# Matches
# 1. Incomplete dates that are missing both the month and day parts (e.g. 2018)
# 2. Other positive integers (e.g. 1, 123, 12345)
RE_YEAR_ONLY = re.compile(r'^0*[1-9][0-9X]*$')


def numeric_date(date, **kwargs):
    """
    Converts the given *date* to a :py:class:`float`.

    Parameters
    ----------
    date
        Date in any of the pre-supported formats or custom format.
        Ambiguous dates must follow the custom format.
    fmt
        Date format string valid for the datetime library
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    min_max_year
        Minimum and maximum year bounds for ambiguous dates.
    ambiguity_resolver
        None: Assume the given date string is an exact date.
        'min': Resolve to minimum of ambiguous range.
        'max': Resolve to maximum of ambiguous range.
        'both': Resolve to (min, max) of ambiguous range.


    Examples
    --------
    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04", fmt="%Y-%m-%d")
    2020.42486...
    >>> import datetime, isodate, treetime
    >>> numeric_date("1W") == treetime.utils.numeric_date(datetime.date.today() - isodate.parse_duration("P1W"))
    True
    """
    try:
        # Put all logic in another function to have the original date string for
        # any errors that may arise.
        return _numeric_date(date, **kwargs)
    except InvalidDateMessage as error:
        raise InvalidDate(date, str(error))


def _numeric_date(date, fmt=None, min_max_year=None, ambiguity_resolver=None):
    # Handle an absolute date as a datetime.date.
    if isinstance(date, datetime.date):
        # Use a treetime utility function to convert the datetime.date to a
        # numeric representation.
        return treetime.utils.numeric_date(date)

    # All other formats are treated as strings.
    date = str(date)

    # Handle an absolute date in numeric format.
    if RE_NUMERIC_DATE.match(date):
        return float(date)

    # Handle an absolute date in a custom date format.
    if fmt:
        # Convert year-only dates to a proper ambiguous format.
        if RE_YEAR_ONLY.match(date):
            date = incomplete_to_ambiguous(fmt, year=date)

        ambiguous_date = AmbiguousDate(date, fmt, min_max_year=min_max_year)

        if ambiguous_date.date_matches_format():
            ambiguous_date.assert_only_less_significant_uncertainty()

            if not 'X' in date:
                # Date is exact and can be parsed given the format.
                date = custom_strptime(date, fmt)
                return treetime.utils.numeric_date(date)

            min_date, max_date = ambiguous_date.range()
            if ambiguity_resolver == 'min':
                return treetime.utils.numeric_date(min_date)
            if ambiguity_resolver == 'max':
                return treetime.utils.numeric_date(max_date)
            if ambiguity_resolver == 'both':
                return (treetime.utils.numeric_date(min_date), treetime.utils.numeric_date(max_date))

    # Finally, handle a backwards-looking relative date in ISO 8601 duration format.
    # Since this is handled last, there is no need to check the format as done for other formats above.
    try:
        date = relative_iso_to_datetime_date(date)
        return treetime.utils.numeric_date(date)
    except (ValueError, isodate.ISO8601Error):
        pass

    # If all formats have been checked and there were no unhandled errors, return None.
    # This happens for "N/A-like" strings (e.g. '', '?') but can also happen if fmt is not properly specified.
    return None


def custom_strptime(date, fmt):
    """datetime.datetime.strptime with custom error messages."""
    try:
        return datetime.datetime.strptime(date, fmt)
    except ValueError as error:
        if "day is out of range for month" in str(error):
            invalid_date_message = str(error)
        elif "does not match format" in str(error):
            invalid_date_message = f"Date is invalid for format {fmt}. Check if any parts are out of bounds."
        else:
            invalid_date_message = f"Unknown error: {error}"
        raise InvalidDateMessage(invalid_date_message) from error


def relative_iso_to_datetime_date(backwards_duration_str: str, from_date: datetime.date = None):
    """Compute datetime.date from a ISO 8601 duration string relative to a specified date.

    Parameters
    ----------
    backwards_duration_str
        ISO 8601 duration string specifying the duration to go back. The 'P' duration designator is optional.
    from_date
        The date to go back from. Default is current date.


    >>> relative_iso_to_datetime_date('5D', from_date=datetime.date(2018, 3, 25))
    datetime.date(2018, 3, 20)
    >>> relative_iso_to_datetime_date('P5D', from_date=datetime.date(2018, 3, 25))
    datetime.date(2018, 3, 20)
    >>> relative_iso_to_datetime_date('5W', from_date=datetime.date(2018, 3, 25))
    datetime.date(2018, 2, 18)
    >>> relative_iso_to_datetime_date('5Y', from_date=datetime.date(2018, 3, 25))
    datetime.date(2013, 3, 25)
    """
    if from_date is None:
        from_date = datetime.date.today()
    if not backwards_duration_str.startswith('P'):
        backwards_duration_str = 'P'+backwards_duration_str
    return from_date - isodate.parse_duration(backwards_duration_str)


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

def get_numerical_dates(metadata:pd.DataFrame, name_col = None, date_col='date', fmt=None, min_max_year=None):
    if not isinstance(metadata, pd.DataFrame):
        raise AugurError("Metadata should be a pandas.DataFrame.")
    if fmt:
        strains = metadata.index.values
        dates = metadata[date_col].apply(
            lambda date: numeric_date(
                date,
                fmt,
                min_max_year,
                ambiguity_resolver='both'
            )
        ).values
    else:
        strains = metadata.index.values
        dates = metadata[date_col].astype(float)
    return dict(zip(strains, dates))

def get_iso_year_week(year, month, day):
    return datetime.date(year, month, day).isocalendar()[:2]


def incomplete_to_ambiguous(fmt, year=None, month=None, day=None):
    """Convert incomplete date parts to an ambiguous date in the given format."""
    if not year:
        year = 'XXXX'
    if not month:
        month = 'XX'
    if not day:
        day = 'XX'

    return (fmt
        .replace("%Y", year)
        .replace("%m", month)
        .replace("%d", day)
    )
