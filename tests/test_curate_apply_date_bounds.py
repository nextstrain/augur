import pytest
from argparse import Namespace
from freezegun import freeze_time
from augur.curate.apply_date_bounds import Record, DataError
from augur.types import DataErrorMethod
from augur.errors import AugurError


def make_args(lower=None, upper=None):
    return Namespace(
        date_field="date",
        lower_bound=lower,
        upper_bound=upper,
        failure_reporting=DataErrorMethod.ERROR_FIRST
    )

@pytest.mark.parametrize(
    "    date                    , lower        , upper        , expected_interval",
    [
        # An exact date is converted to an interval for explicitness.
        ("2020-01-15"            , "2020-01-10" , "2020-01-20" , "2020-01-15/2020-01-15"),

        # A date representing an interval can be bounded.
        ("2020"                  , "2020-02-10" , "2020-05-20" , "2020-02-10/2020-05-20"),
        ("2020-01-01/2020-07-01" , "2020-02-10" , "2020-05-20" , "2020-02-10/2020-05-20"),

        # Bounds can represent intervals too.
        ("2020"                  , "2020-02"    , "2020-05"    , "2020-02-01/2020-05-31"),
    ],
)
def test_date_formats(date, lower, upper, expected_interval):
    """
    Test various date formats in each field.
    """
    record = Record({"date": date, "rootDate": lower, "collectionDate": upper}, 0)
    args = make_args(lower="rootDate", upper="collectionDate")
    assert record.get_bounded_date(args) == expected_interval


@pytest.mark.parametrize(
    "data, lower, upper, expected_interval",
    [
        # An error is shown if it is entirely out of bounds.
        ({"date": "2020"}, "2020-01-10", None, "2020-01-10/2020-12-31"),
        ({"date": "2020"}, None, "2020-01-10", "2020-01-01/2020-01-10"),
    ],
)
def test_constant_bounds(data, lower, upper, expected_interval):
    """
    Test handling of constant bounds.
    """
    record = Record(data, 0)
    args = make_args(lower=lower, upper=upper)
    assert record.get_bounded_date(args) == expected_interval


@pytest.mark.parametrize(
    "data, lower, upper, expected",
    [
        # When both bounds are defined, the date is constructed from the bounds.
        (
            {"date": "XXXX-XX-XX", "rootDate": "2020-01-10", "collectionDate": "2020-01-20"},
            "rootDate",
            "collectionDate",
            "2020-01-10/2020-01-20"
        ),
        # When a single bound is defined, the date is returned unchanged.
        (
            {"date": "XXXX-XX-XX", "collectionDate": "2020-01-20"},
            None,
            "collectionDate",
            "XXXX-XX-XX"
        ),
        (
            {"date": "XXXX-XX-XX", "rootDate": "2020-01-10"},
            "rootDate",
            None,
            "XXXX-XX-XX"
        ),
    ],
)
def test_unknown_date(data, lower, upper, expected):
    """
    Test handling unknown date with both or single bounds.
    """
    record = Record(data, 0)
    args = make_args(lower=lower, upper=upper)
    assert record.get_bounded_date(args) == expected


@pytest.mark.parametrize(
    "data, lower, upper, expected_message_substring",
    [
        ({}, None, None, "Missing date field 'date'"),
        ({"date": "?"}, None, None, "Unable to parse value from 'date' as a date"),
        ({"date": "2020-01-01", "collectionDate": ""}, None, "collectionDate", "Unable to parse value from 'collectionDate' as a date"),

        # An error is shown if it is entirely out of bounds.
        ({"date": "2020-01-05", "rootDate": "2020-01-10"}, "rootDate", None, "earlier than the lower bound"),
        ({"date": "2020-01-15", "collectionDate": "2020-01-10"}, None, "collectionDate", "later than the upper bound"),
    ],
)
def test_data_errors(data, lower, upper, expected_message_substring):
    """
    Test various data errors.
    """
    record = Record(data, 0)
    args = make_args(lower=lower, upper=upper)
    with pytest.raises(DataError) as exc:
        record.get_bounded_date(args)
    assert expected_message_substring in str(exc.value)


@pytest.mark.parametrize(
    "data, lower, upper, expected_message_substring",
    [
        ({"date": "2020-01-15"}, "invalid-bound", None, "Expected --lower-bound to be a field name or date"),
        ({"date": "2020-01-15"}, None, "invalid-bound", "Expected --upper-bound to be a field name or date"),
        ({"date": "2020-01-15", "today": "something"}, None, "today", "'today' is ambiguous as it is both an alias to the current date and a field name"),
    ],
)
def test_user_errors(data, lower, upper, expected_message_substring):
    """
    Test various user errors.
    """
    record = Record(data, 0)
    args = make_args(lower=lower, upper=upper)
    with pytest.raises(AugurError) as exc:
        record.get_bounded_date(args)
    assert expected_message_substring in str(exc.value)


@freeze_time("2020-01-15")
def test_today():
    """
    Test special handing of "today" as an upper bound.
    """
    record = Record({"date": "2020", "rootDate": "2020-01-10"}, 0)
    args = make_args(lower="rootDate", upper="today")
    assert record.get_bounded_date(args) == "2020-01-10/2020-01-15"
