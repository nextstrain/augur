import re

from augur import validate

import pytest


class TestValidate:
    def test_auspice_config_v2_valid(self):
        # Implied assertion that no exceptions are raised
        validate.auspice_config_v2(
            "tests/data/validation/auspice_config_v2__valid.json"
        )

    def test_auspice_config_v2_invalid(self):
        with pytest.raises(
            validate.JsonValidationError, match="'key' is a required property"
        ):
            validate.auspice_config_v2(
                "tests/data/validation/auspice_config_v2__invalid.json"
            )

    def test_export_v2_valid(self):
        # Implied assertion that no exceptions are raised
        validate.export_v2("tests/data/validation/v2_zika__valid.json")

    def test_export_v2_invalid(self):
        try:
            validate.export_v2("tests/data/validation/v2_zika__invalid.json")
        except validate.JsonValidationError as e:
            for i, prop in enumerate(["nuc", "panels", "updated", "version"]):
                assert re.search(f"'{prop}' is a required property", e.errors[i])

    @pytest.mark.parametrize(
        "filename", ["frequencies.json", "entropy.json", "sequences.json"]
    )
    def test_export_v2_bad_filename(self, filename):
        with pytest.raises(
            validate.JsonValidationError,
            match="for the main `augur export v2` JSON only",
        ):
            validate.export_v2(filename)
