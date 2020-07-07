import re

from augur.util_support.metadata_file import MetadataFile

import pytest


@pytest.fixture()
def prepare_file(tmpdir):
    def _prepare_file(contents):
        with open(f"{tmpdir}/metadata.txt", "w") as file:
            file.write(re.sub(r"^\s*", "", contents))

    return _prepare_file


class TestMetadataFile:
    @pytest.mark.parametrize("fname", ["", "/does/not/exist.txt"])
    def test_read_metadata_file_not_found(self, fname):
        with pytest.raises(FileNotFoundError):
            MetadataFile(fname, None).read()

    @pytest.mark.parametrize(
        "query, expected_strains",
        [
            ("location=='colorado'", ["strainA", "strainB"]),
            ("quality=='good' & location=='colorado'", ["strainA"]),
            (None, ["strainA", "strainB", "strainC"]),
        ],
    )
    def test_read_metadata_query(self, tmpdir, prepare_file, query, expected_strains):
        prepare_file(
            """
            strain,location,quality
            strainA,colorado,good
            strainB,colorado,bad
            strainC,nevada,good
            """
        )

        records, columns = MetadataFile(f"{tmpdir}/metadata.txt", query).read()
        assert list(records.keys()) == expected_strains
        assert list(columns) == ["strain", "location", "quality"]

    def test_read_metadata_bad_query(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain,location,quality
            strainA,colorado,good
            strainB,colorado,bad
            strainC,nevada,good
            """
        )

        with pytest.raises(ValueError, match="Error applying pandas query"):
            MetadataFile(f"{tmpdir}/metadata.txt", "age=5").read()

    def test_read_metadata_duplicate_strain(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain,quality
            strainA,good
            strainA,good
            """
        )

        with pytest.raises(ValueError, match="Duplicated strain in metadata: strainA"):
            MetadataFile(f"{tmpdir}/metadata.txt", None).read()

    def test_read_metadata_duplicate_name(self, tmpdir, prepare_file):
        prepare_file(
            """
            name,quality
            nameA,good
            nameA,good
            """
        )

        with pytest.raises(ValueError, match="Duplicated name in metadata: nameA"):
            MetadataFile(f"{tmpdir}/metadata.txt", None).read()

    def test_read_metadata_strain_and_name(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain,name,quality
            strainA,nameA,good
            strainB,nameB,good
            """
        )

        assert MetadataFile(f"{tmpdir}/metadata.txt", None).find_key_type() == "strain"

    def test_read_metadata_no_strain_or_name(self, tmpdir, prepare_file):
        prepare_file(
            """
            location,quality
            colorado,good
            nevada,good
            """
        )

        with pytest.raises(ValueError, match="does not contain"):
            MetadataFile(f"{tmpdir}/metadata.txt", None).read()

    def test_metadata_delimiter_autodetect(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain\tlocation\tquality
            strainA\tcolorado\tgood
            strainB\tnevada\tgood
            """
        )

        records, columns = MetadataFile(f"{tmpdir}/metadata.txt").read()
        assert records == {
            "strainA": {"strain": "strainA", "location": "colorado", "quality": "good"},
            "strainB": {"strain": "strainB", "location": "nevada", "quality": "good"},
        }
        assert list(columns) == ["strain", "location", "quality"]
