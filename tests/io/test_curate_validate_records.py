import pytest
from augur.curate import validate_records
from augur.errors import AugurError


@pytest.fixture
def good_records():
    return [
        {"geo_loc_name": "Canada/Vancouver"},
        {"geo_loc_name": "Canada/Vancouver"},
    ]


@pytest.fixture
def bad_records():
    return [
        {"geo_loc_name": "Canada/Vancouver"},
        {"geo_loc_name2": "Canada/Vancouver"},
    ]


class TestCurateValidateRecords:
    def test_validate_input(self, good_records):
        validated_records = validate_records(good_records, "test_subcmd", True)
        assert list(validated_records) == good_records, "good input records validate"

    def test_validate_output(self, good_records):
        validated_records = validate_records(good_records, "test_subcmd", False)

        assert list(validated_records) == good_records, "good output records validate"

    def test_validate_bad_records(self, bad_records):
        with pytest.raises(AugurError) as e:
            list(validate_records(bad_records, "test_subcmd", True))
        assert str(e.value).startswith(
            "Records do not have the same fields!"
        ), "bad input records throw exception with expected message"

    def test_validate_bad_output(self, bad_records):
        with pytest.raises(AugurError) as e:
            list(validate_records(bad_records, "test_subcmd", False))
        assert str(e.value).startswith(
            "Records do not have the same fields!"
        ), "bad output records throw exception with expected message"
        assert (
          "test_subcmd" in str(e.value)
        ), "bad output records throw exception with subcmd name in the message"
