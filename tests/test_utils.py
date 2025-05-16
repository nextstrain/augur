import json
import numpy as np
from pathlib import Path
from unittest.mock import patch
import pandas as pd

import pytest

from augur import utils
from augur.errors import AugurError


class TestUtils:
    @pytest.mark.parametrize("extension", ["bed","BED"])
    @patch('augur.utils.read_bed_file')
    def test_load_mask_sites_recognizes_bed_file(self, m_read_bed_file, extension):
        """load_mask_sites should handle files that end with .bed with any capitalization as a bed"""
        m_read_bed_file.return_value = [3,4]
        assert utils.load_mask_sites("mask.%s" % extension) == [3,4]
        m_read_bed_file.assert_called_with("mask.%s" % extension)

    @patch('augur.utils.read_mask_file')
    def test_load_mask_sites_recognizes_non_bed_file(self, m_read_mask_file):
        """load_mask_sites should pass any other files to read_mask_file"""
        m_read_mask_file.return_value = [5,6]
        assert utils.load_mask_sites("mask.not_a_bed")
        m_read_mask_file.assert_called_with("mask.not_a_bed")

    def test_read_mask_file_good_input(self, tmpdir):
        """read_mask_file should return a sorted list of unique sites from good mask files"""
        mask_sites = [15,200,34,200,36]
        expected_sites = sorted(set(i - 1 for i in mask_sites))
        mask_file = str(tmpdir / "temp.mask")
        with open(mask_file, "w") as fh:
            fh.write("\n".join(str(i) for i in mask_sites))
        assert utils.read_mask_file(mask_file) == expected_sites

    def test_read_mask_file_bad_lines(self, tmpdir):
        """read_mask_file should fail on bad lines in mask files"""
        mask_file = str(tmpdir / "temp.mask")
        bad_mask_sites = ["1", "#comment", "2"]
        with open(mask_file, "w") as fh:
            fh.write("\n".join(bad_mask_sites))
        with pytest.raises(ValueError):
            utils.read_mask_file(mask_file)

    def test_read_bed_file_good_input(self, tmpdir):
        """read_bed_file should read site ranges from properly formatted bed files"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "SEQ\t2\t6", "SEQ\t3\t4"]
        expected_sites = [2,3,4,5,7]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        assert utils.read_bed_file(bed_file) == expected_sites

    def test_read_bed_file_empty_input(self, tmpdir):
        """read_byd_file should return an empty list given an empty file"""
        bed_file = str(tmpdir / "temp.bed")
        expected_sites = []
        with open(bed_file, "w") as fh:
            fh.write("")
        assert utils.read_bed_file(bed_file) == expected_sites

    def test_read_bed_file_with_header(self, tmpdir):
        """read_bed_file should skip header lines if they exist in bed files"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["CHROM\tSTART\tEND","browser yaddayadda","track moreyadda","SEQ\t7\t8", "SEQ\t2\t5"]
        expected_sites = [2,3,4,7]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        assert utils.read_bed_file(bed_file) == expected_sites

    def test_read_bed_file_with_only_header(self, tmpdir):
        """read_bed_file should skip header lines if they exist in bed files"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["CHROM\tSTART\tEND"]
        expected_sites = []
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        assert utils.read_bed_file(bed_file) == expected_sites

    def test_read_bed_file_with_internal_comment_line(self, tmpdir):
        """read_bed_file should ignore internal comment line"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "# CHROM\tSTART\tEND", "SEQ\t2\t5"]
        expected_sites = [2,3,4,7]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        assert utils.read_bed_file(bed_file) == expected_sites

    def test_read_bed_file_with_internal_header(self, tmpdir):
        """read_bed_file should error out if any other lines are unreadable"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "CHROM\tSTART\tEND", "SEQ\t2\t5"]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        with pytest.raises(AugurError):
            utils.read_bed_file(bed_file)

    def test_read_bed_file_with_mismatched_chrom_values(self, tmpdir):
        """read_bed_file should error out if chrom values don't match"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "SEQ 2\t2\t5"]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        with pytest.raises(AugurError):
            utils.read_bed_file(bed_file)

    def test_read_bed_file_with_bad_data_line(self, tmpdir):
        """read_bed_file should error out if data line has less than 3 values"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "SEQ\t2"]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        with pytest.raises(AugurError):
            utils.read_bed_file(bed_file)

    def test_read_mask_file_drm_file(self, tmpdir):
        """read_mask_file should handle drm files as well"""
        drm_file = str(tmpdir / "temp.drm")
        drm_lines = ["SEQ\t5", "SEQ\t7"]
        expected_sites = [4,6]
        with open(drm_file, "w") as fh:
            fh.write("\n".join(drm_lines))
        assert utils.read_mask_file(drm_file) == expected_sites

    def test_write_json_data_types(self, tmpdir):
        """write_json should be able to serialize various data types."""
        data = {
            'int': np.int64(1),
            'float': np.float64(2.0),
            'array': np.array([3,4,5]),
            'series': pd.Series([6,7,8])
        }
        file = Path(tmpdir) / Path("data.json")
        utils.write_json(data, file, include_version=False)
        with open(file) as f:
            assert json.load(f) == {
                'int': 1,
                'float': 2.0,
                'array': [3,4,5],
                'series': [6,7,8]
            }
