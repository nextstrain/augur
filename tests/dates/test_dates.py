import datetime
import pytest
from freezegun import freeze_time
from augur import dates
from augur.errors import AugurError


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

    def test_get_numerical_date_from_value_not_ambiguous(self):
        assert (dates.get_numerical_date_from_value("2000-03-29", "%Y-%m-%d")
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=3, day=29)), abs=1e-3)
            == pytest.approx(2000.242, abs=1e-3)
        )

    @pytest.mark.parametrize(
        "value",
        [
            "2000-01-XX",
            "2000-01",
            "2000-1",
        ],
    )
    def test_get_numerical_date_from_value_ambiguous_day(self, value):
        min_date, max_date = dates.get_numerical_date_from_value(value, "%Y-%m-%d")
        assert (min_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=1, day=1)), abs=1e-3)
            == pytest.approx(2000.001, abs=1e-3)
        )
        assert (max_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=1, day=31)), abs=1e-3)
            == pytest.approx(2000.083, abs=1e-3)
        )

    @pytest.mark.parametrize(
        "value",
        [
            "2000-001-XX",
            "2000-001",
        ],
    )
    def test_get_numerical_date_from_value_ambiguous_day_error(self, value):
        with pytest.raises(dates.errors.InvalidDate, match="Invalid date '2000-001-XX': Date does not match format `%Y-%m-%d`"):
            dates.get_numerical_date_from_value(value, "%Y-%m-%d")

    @pytest.mark.parametrize(
        "value",
        [
            "2000-01-01",
            "2000-1-1",
            "2000-01-1",
            "2000-1-01",
            "2000-001-01",
            "2000-001-001",
        ],
    )
    def test_get_numerical_date_from_value_iso_like_date(self, value):
        assert (
            dates.get_numerical_date_from_value(value, "%Y-%m-%d")
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=1, day=1)), abs=1e-3)
            == pytest.approx(2000.001, abs=1e-3)
        )

    @pytest.mark.parametrize(
        "value, error",
        [
            ("2000-00-01", "month must be in 1..12"),
            ("2000-13-01", "month must be in 1..12"),
            ("2000-100-01", "month must be in 1..12"),
            ("2000-01-00", "day is out of range for month"),
            ("2000-01-32", "day is out of range for month"),
            ("2000-01-100", "day is out of range for month"),
        ],
    )
    def test_get_numerical_date_from_value_silent_error(self, value, error):
        with pytest.raises(dates.errors.InvalidDate, match=error):
            dates.get_numerical_date_from_value(value, "%Y-%m-%d")

    @pytest.mark.parametrize(
        "value, equivalent",
        [
            ("2000-00", "2000-01"),
            ("2000-13", "2000-12"),
            ("2000-99", "2000-12"),
            ("0000-00", "0001-01"),
        ],
    )
    def test_get_numerical_date_from_value_auto_fix(self, value, equivalent):
        assert dates.get_numerical_date_from_value(value, "%Y-%m-%d") == dates.get_numerical_date_from_value(equivalent, "%Y-%m-%d")

    @freeze_time("2000-02-20")
    @pytest.mark.parametrize(
        "value, equivalent",
        [
            ("2000-03", "2000-02-20"),
        ],
    )
    def test_get_numerical_date_from_value_auto_fix_exact(self, value, equivalent):
        assert dates.get_numerical_date_from_value(value, "%Y-%m-%d") == (dates.get_numerical_date_from_value(equivalent, "%Y-%m-%d"),)*2

    def test_get_numerical_date_from_value_ambiguous_month_and_day(self):
        min_date, max_date = dates.get_numerical_date_from_value("2000-XX-XX", "%Y-%m-%d")
        assert (min_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=1, day=1)), abs=1e-3)
            == pytest.approx(2000.001, abs=1e-3)
        )
        assert (max_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=12, day=31)), abs=1e-3)
            == pytest.approx(2000.999, abs=1e-3)
        )

    @freeze_time("2000-02-20")
    def test_get_numerical_date_from_value_current_day_limit(self):
        min_date, max_date = dates.get_numerical_date_from_value("2000-03-XX", "%Y-%m-%d")
        assert (min_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=2, day=20)), abs=1e-3)
            == pytest.approx(2000.138, abs=1e-3)
        )
        assert (max_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=2, day=20)), abs=1e-3)
            == pytest.approx(2000.138, abs=1e-3)
        )

    def test_get_numerical_date_from_value_range(self):
        assert dates.get_numerical_date_from_value("2019-01-02/2019-03-04", fmt="unused") == (
            pytest.approx(dates.numeric_date(datetime.date(year=2019, month=1, day=2)), abs=1e-3),
            pytest.approx(dates.numeric_date(datetime.date(year=2019, month=3, day=4)), abs=1e-3),
        )

        # Using a numeric date as a bound is not valid.
        assert dates.get_numerical_date_from_value("2019.0/2019-06-01", fmt="unused") == None

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
            dates.get_numerical_dates(metadata, "%Y-%m-%d")
