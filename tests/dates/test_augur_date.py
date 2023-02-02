import pytest
from augur.dates.augur_date import AugurDate
from augur.dates.errors import InvalidDate


class TestExactDate:
    @pytest.mark.parametrize(
        "date_str, expected_str",
        [
            ("2020.0"  , "2020-01-01"),
            ("1.0"     , "0001-01-01"),
            ("0.999"   , "0.999"     ),
            ("0.0"     , "0.0"       ),
            ("-5000.0" , "-5000.0"   ),
            ("-10000.0", "-10000.0"  ),
            ("10000.0" , "10000.0"   ),
        ],
    )
    def test_numeric_date(self, date_str, expected_str):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_str
        assert str(augur_date.max) == expected_str

    @pytest.mark.parametrize(
        "date_str, expected_str",
        [
            ("2020-01-01" , "2020-01-01"),
            ("0001-01-01" , "0001-01-01"),
        ],
    )
    def test_exact_date(self, date_str, expected_str):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_str
        assert str(augur_date.max) == expected_str

    @pytest.mark.parametrize(
        "date_str, expected_min, expected_max",
        [
            ("2000", "2000-01-01", "2000-12-31"),
        ],
    )
    def test_year_only(self, date_str, expected_min, expected_max):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_min
        assert str(augur_date.max) == expected_max

    @pytest.mark.parametrize(
        "date_str, expected_min, expected_max",
        [
            ("2000-01", "2000-01-01", "2000-01-31"),
            ("2000-02", "2000-02-01", "2000-02-29"),
        ],
    )
    def test_year_month(self, date_str, expected_min, expected_max):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_min
        assert str(augur_date.max) == expected_max

    @pytest.mark.parametrize(
        "date_str, expected_min, expected_max",
        [
            ("-1"    , "-1.9999"    , "-1.0"    ),
            ("-5000" , "-5000.9999" , "-5000.0" ),
            ("-10000", "-10000.9999", "-10000.0"),
        ],
    )
    def test_negative_year(self, date_str, expected_min, expected_max):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_min
        assert str(augur_date.max) == expected_max

    @pytest.mark.parametrize(
        "date_str",
        [
            (""),
            ("null"),
            # These cannot be mapped to any default format so they are treated as null.
            ("2000-00"),
            ("2000-13"),
            ("2000-01-00"),
            ("2000-01-32"),
        ],
    )
    def test_null_values(self, date_str):
        augur_date = AugurDate(date_str)
        assert augur_date.is_null

    @pytest.mark.parametrize(
        "date_str, expected_min, expected_max",
        [
            ("[2020-01-01,2020-01-09]", "2020-01-01", "2020-01-09"),
            ("[2020-01-01,2020.10]"   , "2020-01-01", "2020-02-06"),
            ("[-7000.0,-5000.0]"      , "-7000.0"   , "-5000.0"   ),
        ],
    )
    def test_range(self, date_str, expected_min, expected_max):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_min
        assert str(augur_date.max) == expected_max

    @pytest.mark.parametrize(
        "date_str",
        [
            ("[2020-01-09,2020-01-01]"),
            ("[2020-01-09,2020.0]"),
            ("[-5000.0,-7000.0]"),
        ],
    )
    def test_range_order_error(self, date_str):
        with pytest.raises(InvalidDate, match="Range must be in chronological order"):
            AugurDate(date_str)

    @pytest.mark.parametrize(
        "date_str",
        [
            ("[2010,2020]"),
            ("[-7000,-5000]"),
            ("[2020-01-00,2020-01-10]"),
        ],
    )
    def test_range_invalid_error(self, date_str):
        with pytest.raises(InvalidDate, match="Range provided but at least one end is not a valid exact date"):
            AugurDate(date_str)

    @pytest.mark.parametrize(
        "date_str, expected_min, expected_max",
        [
            # Ambiguous month and day
            ("2020"      , "2020-01-01", "2020-12-31"),
            ("2020-XX-XX", "2020-01-01", "2020-12-31"),
            # Ambiguous day
            ("2020-01"   , "2020-01-01", "2020-01-31"),
            ("2020-01-XX", "2020-01-01", "2020-01-31"),
            # Ambiguous day with a range
            ("2020-01-0X", "2020-01-01", "2020-01-09"),
            ("2020-01-X1", "2020-01-01", "2020-01-31"), # TODO: drop support; this isn't really accurate or necessary
            # Ambiguous month with a range
            ("2020-0X"   , "2020-01-01", "2020-09-30"),
            ("2020-0X-XX", "2020-01-01", "2020-09-30"),
            ("2020-1X-XX", "2020-10-01", "2020-12-31"),
            ("2020-X1-XX", "2020-01-01", "2020-12-31"), # TODO: drop support; this isn't really accurate or necessary
            # Ambiguous year with a range
            ("201X"      , "2010-01-01", "2019-12-31"),
        ],
    )
    def test_ambiguous_scalar(self, date_str, expected_min, expected_max):
        augur_date = AugurDate(date_str)
        assert str(augur_date.min) == expected_min
        assert str(augur_date.max) == expected_max

    @pytest.mark.parametrize(
        "date_str, format, expected_min, expected_max",
        [
            # Ambiguous month and day
            ("2020"      , "%Y"      , "2020-01-01", "2020-12-31"),
            ("XX/XX/2020", "%m/%d/%Y", "2020-01-01", "2020-12-31"),
            # Ambiguous day
            ("01/2020"   , "%m/%Y"   , "2020-01-01", "2020-01-31"),
            ("01/XX/2020", "%m/%d/%Y", "2020-01-01", "2020-01-31"),
            # Ambiguous day with a range
            ("01/0X/2020", "%m/%d/%Y", "2020-01-01", "2020-01-09"),
            ("01/X1/2020", "%m/%d/%Y", "2020-01-01", "2020-01-31"),
            # Ambiguous month with a range
            ("0X/2020"   , "%m/%Y"   , "2020-01-01", "2020-09-30"),
            ("0X/XX/2020", "%m/%d/%Y", "2020-01-01", "2020-09-30"),
            ("1X/XX/2020", "%m/%d/%Y", "2020-10-01", "2020-12-31"),
            # Ambiguous year with a range
            ("201X"      , "%Y"      , "2010-01-01", "2019-12-31"),
        ],
    )
    def test_ambiguous_scalar_custom_format(self, date_str, format, expected_min, expected_max):
        augur_date = AugurDate(date_str, possible_formats=[format])
        assert str(augur_date.min) == expected_min
        assert str(augur_date.max) == expected_max
