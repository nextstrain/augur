import pytest
from augur.dates.exact_date import ExactDate
from augur.dates.errors import NotAnExactDate


class TestExactDate:
    @pytest.mark.parametrize(
        "date_str, expected_str",
        [
            ("2020.0" , "2020-01-01"),
            ("1.0"    , "0001-01-01"),
            ("0.999"  , "0.999"     ),
            ("0.0"    , "0.0"       ),
            ("-5000.0", "-5000.0"   ),
        ],
    )
    def test_numeric_dates(self, date_str, expected_str):
        assert str(ExactDate(date_str)) == expected_str

    @pytest.mark.parametrize(
        "date_str, possible_formats, expected_str",
        [
            ("2020-01-01", ["%Y-%m-%d"], "2020-01-01"),
            ("01/01/2020", ["%m/%d/%Y"], "2020-01-01"),
        ],
    )
    def test_custom_formats(self, date_str, possible_formats, expected_str):
        assert str(ExactDate(date_str, possible_formats)) == expected_str

    @pytest.mark.parametrize(
        "date_str, possible_formats, expected_str",
        [
            ("2020-01-01", ["%m/%d/%Y", "%Y-%m-%d"], "2020-01-01"),
            ("01/01/2020", ["%m/%d/%Y", "%Y-%m-%d"], "2020-01-01"),
        ],
    )
    def test_custom_format_order(self, date_str, possible_formats, expected_str):
        assert str(ExactDate(date_str, possible_formats)) == expected_str

    @pytest.mark.parametrize(
        "date_str, possible_formats",
        [
            ("2020-300" , ["%Y-%j"]   ),
            ("2020-43-1", ["%Y-%W-%w"]),
        ],
    )
    def test_unsupported_formats_error(self, date_str, possible_formats):
        """Test date formats that are not yet supported."""
        with pytest.raises(NotAnExactDate):
            ExactDate(date_str, possible_formats)

    @pytest.mark.parametrize(
        "date_str, possible_formats",
        [
            ("2020-01-00", ["%Y-%m-%d"]),
            ("2020-02-30", ["%Y-%m-%d"]),
            ("2020-00-01", ["%Y-%m-%d"]),
            ("2020-13-01", ["%Y-%m-%d"]),
        ],
    )
    def test_out_of_bounds_error(self, date_str, possible_formats):
        """Test dates that are in the correct format but with an out of bounds value."""
        with pytest.raises(NotAnExactDate):
            ExactDate(date_str, possible_formats)

    @pytest.mark.parametrize(
        "date_str, possible_formats",
        [
            ("01/01/2020", ["%Y-%m-%d"]),
            ("2020"      , ["%Y-%m-%d"]),
            ("2020-XX-XX", ["%Y-%m-%d"]),
        ],
    )
    def test_format_mismatch_error(self, date_str, possible_formats):
        with pytest.raises(NotAnExactDate):
            ExactDate(date_str, possible_formats)
