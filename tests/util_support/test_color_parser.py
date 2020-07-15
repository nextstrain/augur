from augur.util_support.color_parser import ColorParser

import pytest


@pytest.fixture()
def prepare_file(tmpdir):
    def _prepare_file(contents):
        with open(f"{tmpdir}/colors.tsv", "w") as file:
            for line in contents.split("\n"):
                file.write(line.strip() + "\n")

    return _prepare_file


class TestColorParser:
    def test_color_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            ColorParser(mapping_filename="/does/not/exist.txt").mapping

    def test_color_file_no_defaults(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain\tstrainA\t#ff0000
            strain\tstrainB\t#00ff00
            strain\tstrainC\t#0000ff
            """
        )

        assert ColorParser(
            mapping_filename=f"{tmpdir}/colors.tsv", use_defaults=False
        ).mapping == {
            "strain": [
                ("straina", "#ff0000"),
                ("strainb", "#00ff00"),
                ("strainc", "#0000ff"),
            ]
        }

    def test_color_file_with_defaults(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain\tstrainA\t#ff0000
            strain\tstrainB\t#00ff00
            strain\tstrainC\t#0000ff
            """
        )

        mapping = ColorParser(mapping_filename=f"{tmpdir}/colors.tsv").mapping
        assert ("africa", "#CEB541") in mapping["region"]
        assert ("north america", "#DC2F24") in mapping["region"]
        assert mapping["strain"] == [
            ("straina", "#ff0000"),
            ("strainb", "#00ff00"),
            ("strainc", "#0000ff"),
        ]

    def test_color_file_bad_entries(self, tmpdir, prepare_file):
        prepare_file(
            """
            strain\tstrainA\t#ff0000
            strain\tstrainB\t#00ff00000000
            strain strainC\t#0000ff
            """
        )

        assert ColorParser(
            mapping_filename=f"{tmpdir}/colors.tsv", use_defaults=False
        ).mapping == {"strain": [("straina", "#ff0000")]}
