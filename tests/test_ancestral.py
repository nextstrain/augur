import pytest
import argparse
import json
from augur import ancestral
from unittest import mock


class TestAncestral:
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
