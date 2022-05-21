import datetime
import pytest
from freezegun import freeze_time
from augur import dates
from augur.dates import InvalidDateFormat, numeric_date
from augur.errors import AugurError


def get_numeric_min_max(date):
    return (
        numeric_date(date, 'min'),
        numeric_date(date, 'max'),
    )


class TestDates:
    @freeze_time("2000-02-20")
    def test_numeric_date(self):
        # Test different representations of February 20, 2000.
        assert dates.numeric_date("2000.138") == pytest.approx(2000.138, abs=1e-3)
        assert dates.numeric_date("2000-02-20") == pytest.approx(2000.138, abs=1e-3)
        assert dates.numeric_date(datetime.date(year=2000, month=2, day=20)) == pytest.approx(2000.138, abs=1e-3)

        # Test relative dates based on freeze_time.
        assert dates.numeric_date("1D") == pytest.approx(2000.135, abs=1e-3)
        assert dates.numeric_date("1W") == pytest.approx(2000.119, abs=1e-3)
        assert dates.numeric_date("1M") == pytest.approx(2000.053, abs=1e-3)
        assert dates.numeric_date("1Y") == pytest.approx(1999.138, abs=1e-3)
        assert dates.numeric_date("1Y1M1W") == pytest.approx(1999.034, abs=1e-3)

    def test_ambiguous_date_to_date_range_not_ambiguous(self):
        assert dates.ambiguous_date_to_date_range("2000-03-29", "%Y-%m-%d") == (
            datetime.date(year=2000, month=3, day=29),
            datetime.date(year=2000, month=3, day=29),
        )

    def test_ambiguous_date_to_date_range_ambiguous_day(self):
        assert dates.ambiguous_date_to_date_range("2000-01-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=1, day=31),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month_and_day(self):
        assert dates.ambiguous_date_to_date_range("2000-XX-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=12, day=31),
        )

    @freeze_time("2000-02-20")
    def test_ambiguous_date_to_date_range_current_day_limit(self):
        assert dates.ambiguous_date_to_date_range("2000-02-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=2, day=1),
            datetime.date(year=2000, month=2, day=20),
        )

    def test_is_date_ambiguous(self):
        """is_date_ambiguous should return true for ambiguous dates and false for valid dates."""
        # Test complete date strings with ambiguous values.
        assert dates.is_date_ambiguous("2019-0X-0X", "any")
        assert dates.is_date_ambiguous("2019-XX-09", "month")
        assert dates.is_date_ambiguous("2019-03-XX", "day")
        assert dates.is_date_ambiguous("201X-03-09", "year")
        assert dates.is_date_ambiguous("20XX-01-09", "month")
        assert dates.is_date_ambiguous("2019-XX-03", "day")
        assert dates.is_date_ambiguous("20XX-01-03", "day")

        # Test incomplete date strings with ambiguous values.
        assert dates.is_date_ambiguous("2019", "any")
        assert dates.is_date_ambiguous("201X", "year")
        assert dates.is_date_ambiguous("2019-XX", "month")
        assert dates.is_date_ambiguous("2019-10", "day")
        assert dates.is_date_ambiguous("2019-XX", "any")
        assert dates.is_date_ambiguous("2019-XX", "day")

        # Test complete date strings without ambiguous dates for the requested field.
        assert not dates.is_date_ambiguous("2019-09-03", "any")
        assert not dates.is_date_ambiguous("2019-03-XX", "month")
        assert not dates.is_date_ambiguous("2019-09-03", "day")
        assert not dates.is_date_ambiguous("2019-XX-XX", "year")

        # Test incomplete date strings without ambiguous dates for the requested fields.
        assert not dates.is_date_ambiguous("2019", "year")
        assert not dates.is_date_ambiguous("2019-10", "month")

    def test_get_numerical_dates_dict_error(self):
        """Using get_numerical_dates with metadata represented as a dict should raise an error."""
        metadata = {
            "example": {
                "strain": "example",
                "date": "2000-03-29"
            }
        }
        with pytest.raises(AugurError):
            dates.get_numerical_dates(metadata)

    def test_ambiguous_day(self):
        """Ambiguous day yields a certain min/max range."""
        date_min, date_max = get_numeric_min_max("2018-01-XX")
        assert date_min == pytest.approx(2018.001, abs=1e-3)
        assert date_max == pytest.approx(2018.083, abs=1e-3)

    def test_missing_day(self):
        """Date without day yields a range equivalent to ambiguous day."""
        date_min, date_max = get_numeric_min_max("2018-01")
        assert date_min == pytest.approx(2018.001, abs=1e-3)
        assert date_max == pytest.approx(2018.083, abs=1e-3)

    def test_ambiguous_month(self):
        """Ambiguous month yields a certain min/max range."""
        date_min, date_max = get_numeric_min_max("2018-XX-XX")
        assert date_min == pytest.approx(2018.001, abs=1e-3)
        assert date_max == pytest.approx(2018.999, abs=1e-3)

    def test_missing_month(self):
        """Date without month/day yields a range equivalent to ambiguous month/day."""
        date_min, date_max = get_numeric_min_max("2018")
        assert date_min == pytest.approx(2018.001, abs=1e-3)
        assert date_max == pytest.approx(2018.999, abs=1e-3)

    def test_numerical_exact_year(self):
        """Numerical year ending in .0 should be interpreted as exact."""
        date_min, date_max = get_numeric_min_max("2018.0")
        assert date_min == pytest.approx(2018.001, abs=1e-3)
        assert date_max == pytest.approx(2018.001, abs=1e-3)

    def test_ambiguous_year(self):
        """Ambiguous year replaces X with 0 (min) and 9 (max)."""
        date_min, date_max = get_numeric_min_max("201X-XX-XX")
        assert date_min == pytest.approx(2010.001, abs=1e-3)
        assert date_max == pytest.approx(2019.999, abs=1e-3)

    def test_ambiguous_year_incomplete_date(self):
        """Ambiguous year without month/day yields a range equivalent to ambiguous month/day counterpart."""
        date_min, date_max = get_numeric_min_max("201X")
        assert date_min == pytest.approx(2010.001, abs=1e-3)
        assert date_max == pytest.approx(2019.999, abs=1e-3)

    @pytest.mark.parametrize(
        "date",
        [
            ("201X-01-01"),
            ("201X-01-XX"),
            ("201X-XX-01"),
            ("2010-XX-01"),
        ],
    )
    def test_invalid_ambiguity_between_date_parts(self, date):
        """Test various forms of invalid ambiguity between date parts."""
        with pytest.raises(InvalidDateFormat):
            get_numeric_min_max(date)

    @pytest.mark.parametrize(
        "date",
        [
            ("20X0"),
            ("20X0-XX-XX"),
            ("2010-X1"),
            ("2010-X1-XX"),
            ("2010-01-X0"),
        ],
    )
    def test_invalid_ambiguity_within_date_part(self, date):
        """Test various forms of invalid ambiguity within date parts."""
        with pytest.raises(InvalidDateFormat):
            get_numeric_min_max(date)

    def test_ambiguous_year_lowercase_x(self):
        """Ambiguous year with a lowercase x raises an error."""
        with pytest.raises(InvalidDateFormat):
            get_numeric_min_max("201x")

    @pytest.mark.parametrize(
        "date",
        [
            ("2018-00-01"),
            ("2018-13-01"),
            ("2018-01-00"),
            ("2018-02-30"),
        ],
    )
    def test_out_of_bounds_error(self, date):
        """Out-of-bounds month/day cannot be parsed."""
        with pytest.raises(InvalidDateFormat):
            get_numeric_min_max(date)

    @pytest.mark.parametrize(
        "date",
        [
            ("-2018-01-01"),
            ("-2018-XX-XX"),
            ("-2018-01"),
            ("-2018"),
        ],
    )
    def test_negative_iso_date_error(self, date):
        """All forms of negative ISO dates are unsupported."""
        with pytest.raises(InvalidDateFormat):
            get_numeric_min_max(date)

    def test_negative_numeric_date(self):
        """Parse negative numeric date."""
        date_min, date_max = get_numeric_min_max("-2018.0")
        assert date_min == pytest.approx(-2018.0, abs=1e-3)
        assert date_max == pytest.approx(-2018.0, abs=1e-3)

    def test_zero_year_incomplete_error(self):
        """Zero year-only date is unsupported."""
        with pytest.raises(InvalidDateFormat):
            get_numeric_min_max("0")

    def test_zero_year_exact(self):
        """Parse the date 0.0."""
        date_min, date_max = get_numeric_min_max("0.0")
        assert date_min == pytest.approx(0.0, abs=1e-3)
        assert date_max == pytest.approx(0.0, abs=1e-3)
