import datetime
from freezegun import freeze_time
import pytest

from augur.utils import ambiguous_date_to_date_range, read_metadata


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

    def test_read_metadata_valid_csv(self, tmpdir):
        valid_metadata = str(tmpdir / "metadata_valid.csv")
        with open(valid_metadata, "w") as fh:
            fh.write("""\
strain,date,region
A/SomeStrain/123/2019,2019-10-01,China
A/AnotherStrain/456/2019,2019-12-01,Germany
""")

        metadata_dict, metadata_columns = read_metadata(valid_metadata)
        assert len(metadata_columns) == 3
        assert "A/SomeStrain/123/2019" in metadata_dict
        assert metadata_dict["A/SomeStrain/123/2019"]["date"] == "2019-10-01"

    def test_read_metadata_valid_tsv(self, tmpdir):
        valid_metadata = str(tmpdir / "metadata_valid.tsv")
        with open(valid_metadata, "w") as fh:
            fh.write("""\
strain	date	region
A/SomeStrain/123/2019	2019-10-01	China
A/AnotherStrain/456/2019	2019-12-01	Germany
""")

        metadata_dict, metadata_columns = read_metadata(valid_metadata)
        assert len(metadata_columns) == 3
        assert "A/SomeStrain/123/2019" in metadata_dict
        assert metadata_dict["A/SomeStrain/123/2019"]["date"] == "2019-10-01"

    @pytest.mark.xfail
    def test_read_metadata_valid_table(self, tmpdir):
        # Tab-delimited content in a file without a ".tsv" extension will not
        # get parsed properly.
        valid_metadata = str(tmpdir / "metadata_valid.tab")
        with open(valid_metadata, "w") as fh:
            fh.write("""\
strain	date	region
A/SomeStrain/123/2019	2019-10-01	China
A/AnotherStrain/456/2019	2019-12-01	Germany
""")

        metadata_dict, metadata_columns = read_metadata(valid_metadata)
        assert len(metadata_columns) == 3
        assert "A/SomeStrain/123/2019" in metadata_dict
        assert metadata_dict["A/SomeStrain/123/2019"]["date"] == "2019-10-01"

    def test_read_metadata_missing_file(self):
        with pytest.raises(FileNotFoundError):
            metadata_dict, metadata_columns = read_metadata("")

    def test_read_metadata_bad_id(self, tmpdir):
        invalid_metadata = str(tmpdir / "metadata_invalid.tsv")
        with open(invalid_metadata, "w") as fh:
            fh.write("""\
invalid_strain_id	date	region
A/SomeStrain/123/2019	2019-10-01	China
A/AnotherStrain/456/2019	2019-12-01	Germany
""")

        with pytest.raises(KeyError):
            metadata_dict, metadata_columns = read_metadata(invalid_metadata)

    def test_read_metadata_with_duplicates(self, tmpdir):
        invalid_metadata = str(tmpdir / "metadata_invalid.tsv")
        with open(invalid_metadata, "w") as fh:
            fh.write("""\
strain	date	region
A/SomeStrain/123/2019	2019-10-01	China
A/SomeStrain/123/2019	2019-10-02	China
A/AnotherStrain/456/2019	2019-12-01	Germany
""")

        with pytest.raises(ValueError):
            metadata_dict, metadata_columns = read_metadata(invalid_metadata)
