import pytest
import augur.io.vcf
from treetime.vcf_utils import read_vcf


@pytest.fixture
def mock_run_shell_command(mocker):
    mocker.patch("augur.io.vcf.run_shell_command")


class TestVCF:
    # The `read_vcf` functionality used to be in an augur module when these
    # tests were originally written but we now use TreeTime's function of the
    # same name. The tests remain here to protect against any unforeseen changes.
    def test_read_vcf_compressed(self):
        seq_keep = list(read_vcf("tests/data/tb_lee_2015.vcf.gz")['sequences'].keys())

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"

    def test_read_vcf_uncompressed(self):
        seq_keep = list(read_vcf("tests/data/tb_lee_2015.vcf")['sequences'].keys())

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"

    def test_write_vcf_compressed_input(self, mock_run_shell_command):
        augur.io.vcf.write_vcf(
            "tests/data/tb_lee_2015.vcf.gz", "output_file.vcf.gz", []
        )

        augur.io.vcf.run_shell_command.assert_called_once_with(
            "vcftools --gzvcf tests/data/tb_lee_2015.vcf.gz --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_input(self, mock_run_shell_command):
        augur.io.vcf.write_vcf(
            "tests/data/tb_lee_2015.vcf", "output_file.vcf.gz", []
        )

        augur.io.vcf.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/data/tb_lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_compressed_output(self, mock_run_shell_command):
        augur.io.vcf.write_vcf(
            "tests/data/tb_lee_2015.vcf", "output_file.vcf.gz", []
        )

        augur.io.vcf.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/data/tb_lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_output(self, mock_run_shell_command):
        augur.io.vcf.write_vcf(
            "tests/data/tb_lee_2015.vcf", "output_file.vcf", []
        )

        augur.io.vcf.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/data/tb_lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )

    def test_write_vcf_dropped_samples(self, mock_run_shell_command):
        augur.io.vcf.write_vcf(
            "tests/data/tb_lee_2015.vcf", "output_file.vcf", ["x", "y", "z"]
        )

        augur.io.vcf.run_shell_command.assert_called_once_with(
            "vcftools --remove-indv x --remove-indv y --remove-indv z --vcf tests/data/tb_lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )
