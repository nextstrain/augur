import pytest
import argparse
import json
import filecmp
from augur import ancestral
from unittest import mock


class TestAncestral:
    """
    Augur's end-to-end tests execute Snakemake pipelines that call a series of augur CLI commands.

    Some of those pipelines call 'augur ancestral', including:
        - add_to_alignment
        - tb
        - tb_drm
        - zika

    This test class reproduces each of those pipelines' calls to 'augur ancestral', using the data that each pipeline passes in.  
    For each build, the input files are stored in tests/data/ancestral/<build name>/in/.

    Each test then compares the file output of that call to a control file that the pipeline has actually generated.  
    For each build, the expected control files are stored in tests/data/ancestral/<build name>/expected/.
    """

    def test_add_to_alignment_run(self, tmpdir):
        """
        Reproduce the behavior of the 'augur ancestral' CLI command, as implemented in the 'add_to_alignment' test build:

        augur ancestral --tree results/tree.nwk --alignment results/aligned.fasta --output results/nt_muts.json --inference joint
        """
        TEST_DATA_DIR = "./tests/data/ancestral/add_to_alignment/"
        TMP_DIR = str(tmpdir + "/")
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/aligned.fasta",
            output_node_data=None,
            output=TMP_DIR + "out/nt_muts.json",
            output_sequences=None,
            inference="joint",
            vcf_reference=None,
            output_vcf=None,
            infer_ambiguous=True,
            keep_overhangs=False,
            version=None,
        )
        result = ancestral.run(args)
        assert result == 0
        with open(TEST_DATA_DIR + "expected/nt_muts.json", "r") as fh:
            nt_muts_expected = json.load(fh)
        with open(TMP_DIR + "out/nt_muts.json", "r") as fh:
            nt_muts_actual = json.load(fh)
        assert nt_muts_actual == nt_muts_expected

    def test_tb_run(self, tmpdir):
        """
        Reproduce the behavior of the 'augur ancestral' CLI command, as implemented in the 'tb' test build:

        augur ancestral --tree results/tree.nwk --alignment results/masked.vcf.gz --output results/nt_muts.json --inference joint --output-vcf results/nt_muts.vcf --vcf-reference data/ref.fasta
        """
        TEST_DATA_DIR = "./tests/data/ancestral/tb/"
        TMP_DIR = str(tmpdir + "/")
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/masked.vcf.gz",
            output_node_data=None,
            output=TMP_DIR + "out/nt_muts.json",
            output_sequences=None,
            inference="joint",
            vcf_reference=TEST_DATA_DIR + "in/ref.fasta",
            output_vcf=TMP_DIR + "out/nt_muts.vcf",
            infer_ambiguous=True,
            keep_overhangs=False,
            version=None,
        )
        result = ancestral.run(args)
        assert result == 0
        with open(TEST_DATA_DIR + "expected/nt_muts.json", "r") as fh:
            nt_muts_json_expected = json.load(fh)
        with open(TMP_DIR + "out/nt_muts.json", "r") as fh:
            nt_muts_json_actual = json.load(fh)
        assert nt_muts_json_actual == nt_muts_json_expected
        # Check the VCF output as well
        assert filecmp.cmp(
            TMP_DIR + "out/nt_muts.vcf",
            TEST_DATA_DIR + "expected/nt_muts.vcf",
            shallow=False,
        )

    def test_tb_drm_run(self, tmpdir):
        """
        Reproduce the behavior of the 'augur ancestral' CLI command, as implemented in the 'tb_drm' test build:

        augur ancestral --tree results/tree.nwk --alignment results/masked.vcf.gz --output results/nt_muts.json --inference joint --output-vcf results/nt_muts.vcf --vcf-reference data/ref.fasta
        """
        TEST_DATA_DIR = "./tests/data/ancestral/tb_drm/"
        TMP_DIR = str(tmpdir + "/")
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/masked.vcf.gz",
            output_node_data=None,
            output=TMP_DIR + "out/nt_muts.json",
            output_sequences=None,
            inference="joint",
            vcf_reference=TEST_DATA_DIR + "in/ref.fasta",
            output_vcf=TMP_DIR + "out/nt_muts.vcf",
            infer_ambiguous=True,
            keep_overhangs=False,
            version=None,
        )
        result = ancestral.run(args)
        assert result == 0
        with open(TEST_DATA_DIR + "expected/nt_muts.json", "r") as fh:
            nt_muts_json_expected = json.load(fh)
        with open(TMP_DIR + "out/nt_muts.json", "r") as fh:
            nt_muts_json_actual = json.load(fh)
        assert nt_muts_json_actual == nt_muts_json_expected
        # Check the VCF output as well
        assert filecmp.cmp(
            TMP_DIR + "out/nt_muts.vcf",
            TEST_DATA_DIR + "expected/nt_muts.vcf",
            shallow=False,
        )

    def test_zika_run(self, tmpdir):
        """
        Reproduce the behavior of the 'augur ancestral' CLI command, as implemented in the 'zika' test build:

        augur ancestral --tree results/tree.nwk --alignment results/aligned.fasta --infer-ambiguous --output results/nt_muts.json --inference joint
        """
        TEST_DATA_DIR = "./tests/data/ancestral/zika/"
        TMP_DIR = str(tmpdir + "/")
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/aligned.fasta",
            output_node_data=None,
            output=TMP_DIR + "out/nt_muts.json",
            output_sequences=None,
            inference="joint",
            vcf_reference=None,
            output_vcf=None,
            infer_ambiguous=True,
            keep_overhangs=False,
            version=None,
        )
        result = ancestral.run(args)
        assert result == 0
        with open(TEST_DATA_DIR + "expected/nt_muts.json", "r") as fh:
            nt_muts_expected = json.load(fh)
        with open(TMP_DIR + "out/nt_muts.json", "r") as fh:
            nt_muts_actual = json.load(fh)
        assert nt_muts_actual == nt_muts_expected
