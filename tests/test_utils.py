import datetime

from augur.utils import ambiguous_date_to_date_range

from freezegun import freeze_time


class TestUtils:
    def test_ambiguous_date_to_date_range_not_ambiguous(self):
        assert ambiguous_date_to_date_range("2000-03-29", "%Y-%m-%d") == (
            datetime.date(year=2000, month=3, day=29),
            datetime.date(year=2000, month=3, day=29),
        )

    def test_ambiguous_date_to_date_range_ambiguous_day(self):
        assert ambiguous_date_to_date_range("2000-01-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=1, day=31),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month(self):
        assert ambiguous_date_to_date_range("2000-XX-5", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=5),
            datetime.date(year=2000, month=12, day=5),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month_and_day(self):
        assert ambiguous_date_to_date_range("2000-XX-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=12, day=31),
        )

    @freeze_time("2000-02-20")
    def test_ambiguous_date_to_date_range_current_day_limit(self):
        assert ambiguous_date_to_date_range("2000-02-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=2, day=1),
            datetime.date(year=2000, month=2, day=20),
        )
