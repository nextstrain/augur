import datetime
from unittest.mock import patch

import pytest
from freezegun import freeze_time

from augur import utils

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

    def test_read_bed_file_with_header(self, tmpdir):
        """read_bed_file should skip header lines if they exist in bed files"""
        bed_file = str(tmpdir / "temp.bed")
        bed_lines = ["CHROM\tSTART\tEND","SEQ\t7\t8", "SEQ\t2\t5"]
        expected_sites = [2,3,4,7]
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
    
    def test_read_mask_file_drm_file(self, tmpdir):
        """read_mask_file should handle drm files as well"""
        drm_file = str(tmpdir / "temp.drm")
        drm_lines = ["SEQ\t5", "SEQ\t7"]
        expected_sites = [4,6]
        with open(drm_file, "w") as fh:
            fh.write("\n".join(drm_lines))
        assert utils.read_mask_file(drm_file) == expected_sites

    def test_read_metadata_no_filename(self):
        expected_result = {}, []
        assert utils.read_metadata('') == expected_result

    def test_read_metadata_bad_filename(self, tmpdir):
        bad_filename = "incorrect_path.tsv"
        expected_result = {}, []
        
        meta_fn = str(tmpdir / "metadata.tsv")
        meta_lines = ["strain\tlocation\tquality",
                      "c_good\tcolorado\tgood",
                      "n_bad\tnevada\tbad",]
        with open(meta_fn, "w") as fh:
            fh.write("\n".join(meta_lines))
        
        utils.read_metadata(bad_filename, query='quality=="good" & location=="colorado"')
        assert utils.read_metadata('') == expected_result

    def test_read_metadata_with_good_query(self, tmpdir):
        meta_fn = str(tmpdir / "metadata.tsv")
        meta_lines = ["strain\tlocation\tquality",
                      "c_good\tcolorado\tgood",
                      "c_bad\tcolorado\tbad",
                      "n_good\tnevada\tgood"]
        with open(meta_fn, "w") as fh:
            fh.write("\n".join(meta_lines))
        meta_dict, _ = utils.read_metadata(meta_fn, query='quality=="good" & location=="colorado"')
        assert len(meta_dict) == 1
        assert "c_good" in meta_dict

    def test_read_metadata_bad_query(self, tmpdir):
        meta_fn = str(tmpdir / "metadata.tsv")
        meta_lines = ["strain\tlocation\tquality",
                      "c_good\tcolorado\tgood",
                      "n_bad\tnevada\tbad",]
        with open(meta_fn, "w") as fh:
            fh.write("\n".join(meta_lines))
        with pytest.raises(SystemExit):
            utils.read_metadata(meta_fn, query='badcol=="goodval"')    
    def test_read_metadata_duplicate_strain(self, tmpdir):
        meta_fn = str(tmpdir / "metadata.tsv")
        meta_lines = ["strain\tlocation\tquality",
                      "c_good\tcolorado\tgood",
                      "c_bad\tcolorado\tbad",
                      "c_bad\tcolorado\tbad",
                      "n_good\tnevada\tgood"]
        with open(meta_fn, "w") as fh:
            fh.write("\n".join(meta_lines))
        with pytest.raises(ValueError): 
            utils.read_metadata(meta_fn, query='location=="colorado"')

    def test_read_metadata_duplicate_name(self, tmpdir):
        meta_fn = str(tmpdir / "metadata.tsv")
        meta_lines = ["name\tlocation\tquality",
                      "c_good\tcolorado\tgood",
                      "c_bad\tcolorado\tbad",
                      "c_bad\tcolorado\tbad",
                      "n_good\tnevada\tgood"]
        with open(meta_fn, "w") as fh:
            fh.write("\n".join(meta_lines))
        
        with pytest.raises(ValueError): 
            utils.read_metadata(meta_fn)
        
    def test_read_metadata_no_strain_or_name(self, tmpdir):
        expected_result = {}, []
        
        meta_fn = str(tmpdir / "metadata.tsv")
        meta_lines = ["location\tquality",
                      "colorado\tgood",
                      "colorado\tbad",
                      "colorado\tbad",
                      "nevada\tgood"]
        with open(meta_fn, "w") as fh:
            fh.write("\n".join(meta_lines))
        
        assert utils.read_metadata('') == expected_result
