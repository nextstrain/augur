import pytest
import argparse
import json
import filecmp
from augur import ancestral
from unittest import mock


class TestAncestral:
    def test_add_to_alignment_run(self):
        TEST_DATA_DIR = "./tests/data/ancestral/add_to_alignment/"
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/aligned.fasta",
            output_node_data=None,
            output=TEST_DATA_DIR + "out/nt_muts.json",
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
        nt_muts_expected = json.load(open(TEST_DATA_DIR + "expected/nt_muts.json", "r"))
        nt_muts_actual = json.load(open(TEST_DATA_DIR + "out/nt_muts.json", "r"))
        assert nt_muts_actual == nt_muts_expected

    def test_tb_run(self):
        TEST_DATA_DIR = "./tests/data/ancestral/tb/"
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/masked.vcf.gz",
            output_node_data=None,
            output=TEST_DATA_DIR + "out/nt_muts.json",
            output_sequences=None,
            inference="joint",
            vcf_reference=TEST_DATA_DIR + "in/ref.fasta",
            output_vcf=TEST_DATA_DIR + "out/nt_muts.vcf",
            infer_ambiguous=True,
            keep_overhangs=False,
            version=None,
        )
        result = ancestral.run(args)
        assert result == 0
        nt_muts_json_expected = json.load(
            open(TEST_DATA_DIR + "expected/nt_muts.json", "r")
        )
        nt_muts_json_actual = json.load(open(TEST_DATA_DIR + "out/nt_muts.json", "r"))
        assert nt_muts_json_actual == nt_muts_json_expected
        # Check the VCF output as well
        assert filecmp.cmp(
            TEST_DATA_DIR + "out/nt_muts.vcf",
            TEST_DATA_DIR + "expected/nt_muts.vcf",
            shallow=False,
        )

    def test_tb_drm_run(self):
        TEST_DATA_DIR = "./tests/data/ancestral/tb_drm/"
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/masked.vcf.gz",
            output_node_data=None,
            output=TEST_DATA_DIR + "out/nt_muts.json",
            output_sequences=None,
            inference="joint",
            vcf_reference=TEST_DATA_DIR + "in/ref.fasta",
            output_vcf=TEST_DATA_DIR + "out/nt_muts.vcf",
            infer_ambiguous=True,
            keep_overhangs=False,
            version=None,
        )
        result = ancestral.run(args)
        assert result == 0
        nt_muts_json_expected = json.load(
            open(TEST_DATA_DIR + "expected/nt_muts.json", "r")
        )
        nt_muts_json_actual = json.load(open(TEST_DATA_DIR + "out/nt_muts.json", "r"))
        assert nt_muts_json_actual == nt_muts_json_expected
        # Check the VCF output as well
        assert filecmp.cmp(
            TEST_DATA_DIR + "out/nt_muts.vcf",
            TEST_DATA_DIR + "expected/nt_muts.vcf",
            shallow=False,
        )

    def test_zika_run(self):
        TEST_DATA_DIR = "./tests/data/ancestral/zika/"
        args = argparse.Namespace(
            tree=TEST_DATA_DIR + "in/tree.nwk",
            alignment=TEST_DATA_DIR + "in/aligned.fasta",
            output_node_data=None,
            output=TEST_DATA_DIR + "out/nt_muts.json",
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
        nt_muts_expected = json.load(open(TEST_DATA_DIR + "expected/nt_muts.json", "r"))
        nt_muts_actual = json.load(open(TEST_DATA_DIR + "out/nt_muts.json", "r"))
        assert nt_muts_actual == nt_muts_expected
