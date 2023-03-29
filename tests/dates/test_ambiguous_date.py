import datetime

from augur.dates import ambiguous_date
from augur.dates.ambiguous_date import AmbiguousDate
from augur.dates.errors import InvalidDate

from freezegun import freeze_time
import pytest


class TestAmbiguousDate:
    @freeze_time("2111-05-05")
    @pytest.mark.parametrize(
        "date_str, expected_range",
        [
            ("2000-01-01", (datetime.date(2000, 1, 1), datetime.date(2000, 1, 1))),
            ("2000-02-XX", (datetime.date(2000, 2, 1), datetime.date(2000, 2, 29))),
            ("2000-XX-XX", (datetime.date(2000, 1, 1), datetime.date(2000, 12, 31))),
        ],
    )
    def test_range(self, date_str, expected_range):
        assert AmbiguousDate(date_str).range() == expected_range

    @pytest.mark.parametrize(
        "date_str, fmt",
        [
            ("2005-02-XX", "%Y-%m-%d"),
            ("2005/02/XX", "%Y/%m/%d"),
            ("2005-XX-02", "%Y-%d-%m"),
            ("200502XX", "%Y%m%d"),
        ],
    )
    def test_range_separators(self, date_str, fmt):
        assert AmbiguousDate(date_str, fmt=fmt).range() == (
            datetime.date(2005, 2, 1),
            datetime.date(2005, 2, 28),
        )

    @pytest.mark.parametrize(
        "date_str, expected_components",
        [
            ("2000-01-01", {"Y": "2000", "m": "01", "d": "01"}),
            ("2000-01-XX", {"Y": "2000", "m": "01", "d": "XX"}),
            ("2000-XX-XX", {"Y": "2000", "m": "XX", "d": "XX"}),
        ],
    )
    def test_uncertain_date_components(self, date_str, expected_components):
        assert (
            AmbiguousDate(date_str).uncertain_date_components == expected_components
        )

    def test_uncertain_date_components_error(self):
        with pytest.raises(InvalidDate, match="Date does not match format"):
            AmbiguousDate("5-5-5-5-5").uncertain_date_components

    @pytest.mark.parametrize(
        "date_str, min_or_max, expected",
        [
            ("2000", "min", 2000),
            ("2000", "max", 2000),
            ("200X", "min", 2000),
            ("200X", "max", 2009),
            ("20X0", "max", 2090),
            ("X000", "max", 9000),
            ("XXXX", "min", 1),
            ("XXXX", "max", 9999),
        ],
    )
    def test_resolve_uncertain_int(self, date_str, min_or_max, expected):
        assert (
            ambiguous_date.resolve_uncertain_int(date_str, min_or_max) == expected
        )

    @pytest.mark.parametrize(
        "date_str, expected_error",
        [
            ("200X-01-01", "so month and day must also be uncertain"),
            ("2000-XX-01", "so day must also be uncertain"),
        ],
    )
    def test_assert_only_less_significant_uncertainty(self, date_str, expected_error):
        with pytest.raises(InvalidDate, match=expected_error):
            AmbiguousDate(date_str)

    @freeze_time("2020-05-05")
    @pytest.mark.parametrize(
        "date_str, min_max_year, expected_range",
        [
            ("20XX-XX-XX", [2000, 2010], (datetime.date(2000, 1, 1), datetime.date(2010, 12, 31))),

            # This option does not apply when the year is exact. However, a separate limit of the current date is applied.
            ("2020-XX-XX", [2000, 2010], (datetime.date(2020, 1, 1), datetime.date(2020, 5, 5))),
            ("2020-12-XX", [2000, 2010], (datetime.date(2020, 5, 5), datetime.date(2020, 5, 5))),
            ("2020-12-01", [2000, 2010], (datetime.date(2020, 5, 5), datetime.date(2020, 5, 5))),

            # The upper bound is in the future, which is valid. However, a separate limit of the current date is applied.
            ("20XX-XX-XX", [2010, 2030], (datetime.date(2010, 1, 1), datetime.date(2020, 5, 5))),

            # When there is no upper bound, a date can appear in the future. However, a separate limit of the current date is applied.
            ("20XX-XX-XX", [2010], (datetime.date(2010, 1, 1), datetime.date(2020, 5, 5))),
        ],
    )
    def test_min_max_year(self, date_str, min_max_year, expected_range):
        assert (
            AmbiguousDate(date_str).range(min_max_year=min_max_year) == expected_range
        )

    @freeze_time("2020-05-05")
    @pytest.mark.parametrize(
        "date_str, min_max_year",
        [
            ("20XX-XX-XX", [1950, 1960]),
            ("19XX-XX-XX", [2000, 2010]),
        ],
    )
    def test_min_max_year_date_error(self, date_str, min_max_year):
        with pytest.raises(InvalidDate, match="Not possible for date to fall within bounds"):
            AmbiguousDate(date_str).range(min_max_year=min_max_year)
