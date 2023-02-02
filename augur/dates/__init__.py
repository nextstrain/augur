import argparse
import datetime
import isodate
import pandas as pd
import treetime.utils
from augur.errors import AugurError

from .augur_date import AugurDate
from .errors import InvalidDate

SUPPORTED_DATE_HELP_TEXT = None


def numeric_date(date, possible_formats=None, ambiguity_resolver=None):
    # date is a datetime.date
    if isinstance(date, datetime.date):
        return treetime.utils.numeric_date(date)

    # If parse-able as a duration, return the relative date.
    try:
        return treetime.utils.numeric_date(relative_iso_to_datetime_date(date))
    except (ValueError, isodate.ISO8601Error):
        pass

    # Handle custom formats.
    augur_date = AugurDate(date, possible_formats=possible_formats)
    min_date, max_date = augur_date.range()
    if ambiguity_resolver == None:
        if min_date != max_date:
            raise AugurError("Date is ambiguous but no method to resolve ambiguity was provided.")
        return min_date.as_numeric
    if ambiguity_resolver == 'min':
        return min_date.as_numeric
    if ambiguity_resolver == 'max':
        return max_date.as_numeric
    if ambiguity_resolver == 'both':
        return (min_date.as_numeric, max_date.as_numeric)


def relative_iso_to_datetime_date(backwards_duration_str: str, from_date: datetime.date = None):
    """Compute datetime.date from a ISO 8601 duration string relative to a specified date.
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
    return from_date - isodate.parse_duration(backwards_duration_str)


def numeric_date_type(date):
    """Wraps numeric_date() for argparse usage.

    This raises an ArgumentTypeError from InvalidDateFormat exceptions, otherwise the custom exception message won't be shown in console output due to:
    https://github.com/python/cpython/blob/5c4d1f6e0e192653560ae2941a6677fbf4fbd1f2/Lib/argparse.py#L2503-L2513
    """
    try:
        return numeric_date(date, ambiguity_resolver="min")
    except InvalidDate as error:
        raise argparse.ArgumentTypeError(str(error)) from error


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
    date = AugurDate(date)
    if date.is_null:
        return True
    date_min = date.min.as_datetime
    date_max = date.max.as_datetime

    # Determine ambiguity hierarchically such that, for example, an ambiguous
    # month implicates an ambiguous day even when day information is available.
    # TODO: remove the above comment with an explanation that it isn't appropriate when using AmbiguousDate
    if ambiguous_by in {"any", "day"}:
        return date_min != date_max
    if ambiguous_by in {"any", "month", "day"}:
        return date_min.month != date_max.month
    if ambiguous_by in {"any", "day", "month", "year"}:
        return date_min.year != date_max.year


def get_numerical_date_from_value(value, fmt=None, min_max_year=None):
    # TODO: min_max_year
    possible_formats = None
    if fmt:
        possible_formats = [fmt]
    date = AugurDate(value, possible_formats=possible_formats)
    if date.is_null:
        return None
    if date.min == date.max:
        return date.min.as_numeric
    else:
        return (date.min.as_numeric, date.max.as_numeric)


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
    pass
