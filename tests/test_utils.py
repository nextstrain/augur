from pathlib import Path
from unittest.mock import patch

import pytest

from augur import utils
from test_filter_pandas import write_metadata


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

    def test_read_strains(self, tmpdir):
        # Write one list of filenames with some unnecessary whitespace.
        strains1 = Path(tmpdir) / Path("strains1.txt")
        with open(strains1, "w") as oh:
            oh.write("strain1 # this is an inline comment about strain 1\nstrain2\n   # this is a comment preceded by whitespace.\n")

        # Write another list of filenames with a comment.
        strains2 = Path(tmpdir) / Path("strains2.txt")
        with open(strains2, "w") as oh:
            oh.write("# this is a comment. ignore this.\nstrain2\nstrain3\n")

        strains = utils.read_strains(strains1, strains2)
        assert len(strains) == 3
        assert "strain1" in strains

    def test_read_metadata(self, tmpdir):
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_2","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        utils.read_metadata(meta_fn, as_data_frame=True)
        # duplicates SEQ_1 raises ValueError
        meta_fn = write_metadata(tmpdir, (("strain","location","quality"),
                                          ("SEQ_1","colorado","good"),
                                          ("SEQ_1","colorado","bad"),
                                          ("SEQ_3","nevada","good")))
        with pytest.raises(ValueError) as e_info:
            utils.read_metadata(meta_fn, as_data_frame=True)
        assert str(e_info.value) == "Duplicated strain in metadata: SEQ_1"
