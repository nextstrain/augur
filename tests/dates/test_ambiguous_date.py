import datetime

from augur.dates import ambiguous_date
from augur.dates.ambiguous_date import AmbiguousDate
from augur.dates.errors import InvalidDate

from freezegun import freeze_time
import pytest

from augur.dates.exact_date import ExactDate


class TestAmbiguousDate:
    @freeze_time("2111-05-05")
    @pytest.mark.parametrize(
        "date_str, expected_range",
        [
            ("2000-01-01", ("2000-01-01", "2000-01-01")),
            ("2000-02-XX", ("2000-02-01", "2000-02-29")),
            ("2000-XX-XX", ("2000-01-01", "2000-12-31")),
        ],
    )
    def test_range(self, date_str, expected_range):
        format = '%Y-%m-%d'
        assert AmbiguousDate(date_str, fmt=format).range() == (
            ExactDate(expected_range[0], [format]),
            ExactDate(expected_range[1], [format]),
        )

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
            ExactDate('2005-02-01', ["%Y-%m-%d"]),
            ExactDate('2005-02-28', ["%Y-%m-%d"]),
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
