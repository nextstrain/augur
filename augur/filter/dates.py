from functools import lru_cache
import treetime.utils
from augur.dates import get_numerical_date_from_value
from augur.dates.errors import InvalidDate
from augur.errors import AugurError
from augur.filter.debug import add_debugging
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.sqlite3 import Sqlite3Database, sanitize_identifier
from . import constants


@add_debugging
def parse_dates():
    """Validate dates and create a date table."""
    # First, determine if there is a date column.
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_has_date_column = METADATA_DATE_COLUMN in db.columns(constants.METADATA_TABLE)

    if metadata_has_date_column:
        # Check dates for errors in Python since error handling is non-trivial
        # with SQLite3 user-defined functions.
        _validate_dates()

        _create_date_table_from_metadata()
    else:
        # Create a placeholder table so later JOINs to this table will not break.
        _create_empty_date_table()

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.create_primary_index(constants.DATE_TABLE, constants.ID_COLUMN)


def _validate_dates():
    """Query metadata for dates and error upon any invalid dates."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT {METADATA_DATE_COLUMN}
            FROM {constants.METADATA_TABLE}
        """)
        for row in result:
            try:
                get_numerical_date_from_value(str(row[METADATA_DATE_COLUMN]), fmt='%Y-%m-%d')
            except InvalidDate as error:
                raise AugurError(error)


def _create_date_table_from_metadata():
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        # Register SQLite3 user-defined functions.
        db.connection.create_function(get_year.__name__ , 1, get_year)
        db.connection.create_function(get_month.__name__, 1, get_month)
        db.connection.create_function(get_day.__name__  , 1, get_day)
        db.connection.create_function(get_week.__name__ , 1, get_week)
        db.connection.create_function(try_get_numeric_date_min.__name__, 1, try_get_numeric_date_min)
        db.connection.create_function(try_get_numeric_date_max.__name__, 1, try_get_numeric_date_max)

        db.connection.execute(f"""CREATE TABLE {constants.DATE_TABLE} AS
            SELECT
                {sanitize_identifier(metadata_id_column)} AS {constants.ID_COLUMN},
                {get_year.__name__}({METADATA_DATE_COLUMN})  AS {constants.DATE_YEAR_COLUMN},
                {get_month.__name__}({METADATA_DATE_COLUMN}) AS {constants.DATE_MONTH_COLUMN},
                {get_day.__name__}({METADATA_DATE_COLUMN})   AS {constants.DATE_DAY_COLUMN},
                {get_week.__name__}({METADATA_DATE_COLUMN})  AS {constants.DATE_WEEK_COLUMN},
                {try_get_numeric_date_min.__name__}({METADATA_DATE_COLUMN}) AS {constants.NUMERIC_DATE_MIN_COLUMN},
                {try_get_numeric_date_max.__name__}({METADATA_DATE_COLUMN}) AS {constants.NUMERIC_DATE_MAX_COLUMN}
            FROM {constants.METADATA_TABLE}
        """)

        # Remove user-defined functions.
        db.connection.create_function(get_year.__name__ , 1, None)
        db.connection.create_function(get_month.__name__, 1, None)
        db.connection.create_function(get_day.__name__  , 1, None)
        db.connection.create_function(get_week.__name__ , 1, None)
        db.connection.create_function(try_get_numeric_date_min.__name__, 1, None)
        db.connection.create_function(try_get_numeric_date_max.__name__, 1, None)


def _create_empty_date_table():
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        metadata_id_column = db.get_primary_index(constants.METADATA_TABLE)

    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.connection.execute(f"""CREATE TABLE {constants.DATE_TABLE} AS
        SELECT
            {sanitize_identifier(metadata_id_column)} AS {constants.ID_COLUMN},
            NULL AS {METADATA_DATE_COLUMN},
            NULL AS {constants.DATE_YEAR_COLUMN},
            NULL AS {constants.DATE_MONTH_COLUMN},
            NULL AS {constants.DATE_DAY_COLUMN},
            NULL AS {constants.DATE_WEEK_COLUMN},
            NULL AS {constants.NUMERIC_DATE_MIN_COLUMN},
            NULL AS {constants.NUMERIC_DATE_MAX_COLUMN}
        FROM {constants.METADATA_TABLE}
    """)


CACHE_SIZE = 8192
# Some functions below use @lru_cache to minimize redundant operations on large
# datasets that are likely to have multiple entries with the same date value.


@lru_cache(maxsize=CACHE_SIZE)
def get_year(date):
    """Get the year from a date.
    This function is intended to be registered as a user-defined function in sqlite3. As such, it will not raise any errors.
    """
    try:
        date_min, date_max = datetime_range(date)
    except:
        return None

    if date_min.year == date_max.year:
        return f"{date_min.year}"
    return None


@lru_cache(maxsize=CACHE_SIZE)
def get_month(date):
    """Get the month from a date.
    This function is intended to be registered as a user-defined function in sqlite3. As such, it will not raise any errors.
    """
    try:
        date_min, date_max = datetime_range(date)
    except:
        return None

    if date_min.year == date_max.year and date_min.month == date_max.month:
        return f"{date_min.year}-{str(date_min.month).zfill(2)}"
    return None


# FIXME: remove this?
@lru_cache(maxsize=CACHE_SIZE)
def get_day(date):
    """Get the day from a date.
    This function is intended to be registered as a user-defined function in sqlite3. As such, it will not raise any errors.
    """
    try:
        date_min, date_max = datetime_range(date)
    except:
        return None

    if date_min == date_max:
        return f"{date_min.day}"
    return None


@lru_cache(maxsize=CACHE_SIZE)
def get_week(date):
    """Get the year and week from a date.
    This function is intended to be registered as a user-defined function in sqlite3. As such, it will not raise any errors.
    """
    try:
        date_min, date_max = datetime_range(date)
    except:
        return None

    if date_min == date_max:
        year, week = date_min.isocalendar()[:2]
        return f"{year}-{str(week).zfill(2)}"
    return None


def datetime_range(date):
    numeric_min = get_numerical_date_from_value(date, fmt='%Y-%m-%d', ambiguity_resolver='min')
    numeric_max = get_numerical_date_from_value(date, fmt='%Y-%m-%d', ambiguity_resolver='max')
    date_min = treetime.utils.datetime_from_numeric(numeric_min)
    date_max = treetime.utils.datetime_from_numeric(numeric_max)
    return (date_min, date_max)


@lru_cache(maxsize=CACHE_SIZE)
def try_get_numeric_date_min(date):
    """Get the numeric date from any supported date format, taking the minimum possible value if ambiguous.
    This function is intended to be registered as a user-defined function in sqlite3. As such, it will not raise any errors.
    """
    try:
        return get_numerical_date_from_value(date, fmt='%Y-%m-%d', ambiguity_resolver='min')
    except:
        return None


@lru_cache(maxsize=CACHE_SIZE)
def try_get_numeric_date_max(date):
    """Get the numeric date from any supported date format, taking the maximum possible value if ambiguous.
    This function is intended to be registered as a user-defined function in sqlite3. As such, it will not raise any errors.
    """
    try:
        return get_numerical_date_from_value(date, fmt='%Y-%m-%d', ambiguity_resolver='max')
    except:
        return None
