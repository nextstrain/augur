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

    def test_get_numerical_date_from_value_ambiguous_day(self):
        min_date, max_date = dates.get_numerical_date_from_value("2000-01-XX", "%Y-%m-%d")
        assert (min_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=1, day=1)), abs=1e-3)
            == pytest.approx(2000.001, abs=1e-3)
        )
        assert (max_date
            == pytest.approx(dates.numeric_date(datetime.date(year=2000, month=1, day=31)), abs=1e-3)
            == pytest.approx(2000.083, abs=1e-3)
        )

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
