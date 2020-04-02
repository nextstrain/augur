import datetime
import pytest

from augur import utils

from freezegun import freeze_time
from unittest.mock import patch


class TestUtils:
    def test_ambiguous_date_to_date_range_not_ambiguous(self):
        assert utils.ambiguous_date_to_date_range("2000-03-29", "%Y-%m-%d") == (
            datetime.date(year=2000, month=3, day=29),
            datetime.date(year=2000, month=3, day=29),
        )

    def test_ambiguous_date_to_date_range_ambiguous_day(self):
        assert utils.ambiguous_date_to_date_range("2000-01-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=1, day=31),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month(self):
        assert utils.ambiguous_date_to_date_range("2000-XX-5", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=5),
            datetime.date(year=2000, month=12, day=5),
        )

    def test_ambiguous_date_to_date_range_ambiguous_month_and_day(self):
        assert utils.ambiguous_date_to_date_range("2000-XX-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=1, day=1),
            datetime.date(year=2000, month=12, day=31),
        )

    @freeze_time("2000-02-20")
    def test_ambiguous_date_to_date_range_current_day_limit(self):
        assert utils.ambiguous_date_to_date_range("2000-02-XX", "%Y-%m-%d") == (
            datetime.date(year=2000, month=2, day=1),
            datetime.date(year=2000, month=2, day=20),
        )

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
        """read_mask_file should ignore bad lines in mask files"""
        mask_file = str(tmpdir / "temp.mask")
        bad_mask_sites = ["1", "#comment", "2"]
        expected_sites = [0,1]
        with open(mask_file, "w") as fh:
            fh.write("\n".join(bad_mask_sites))
        assert utils.read_mask_file(mask_file) == expected_sites
    
    def test_read_bed_file_good_input(self, tmpdir):
        """read_bed_file should read site ranges from properly formatted bed files"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "SEQ\t2\t5", "SEQ\t3\t4"]
        expected_sites = [2,3,4,5,7,8]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        assert utils.read_bed_file(bed_file) == expected_sites
    
    def test_read_bed_file_with_header(self, tmpdir):
        """read_bed_file should skip header lines if they exist in bed files"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["CHROM\tSTART\tEND","SEQ\t7\t8", "SEQ\t2\t5"]
        expected_sites = [2,3,4,5,7,8]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        assert utils.read_bed_file(bed_file) == expected_sites
    
    def test_read_bed_file_with_bad_lines(self, tmpdir):
        """read_bed_file should error out if any other lines are unreadable"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["SEQ\t7\t8", "CHROM\tSTART\tEND", "SEQ\t2\t5"]
        with open(bed_file, "w") as fh:
            fh.write("\n".join(bed_lines))
        with pytest.raises(Exception):
            utils.read_bed_file(bed_file)
