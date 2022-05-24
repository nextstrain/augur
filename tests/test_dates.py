import datetime
from freezegun import freeze_time
from augur import dates


class TestDates:
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
