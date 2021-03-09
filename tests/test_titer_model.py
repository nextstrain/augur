"""Tests for titer model module.
"""
import pytest

from augur.reconstruct_sequences import load_alignments
from augur.titer_model import SubstitutionModel, TreeModel
from augur.utils import read_tree


@pytest.fixture
def tree():
    """Returns an annotated Bio.Phylo tree.
    """
    tree = read_tree("tests/data/titer_model/tree_raw.nwk")
    return tree


@pytest.fixture
def alignments():
    alignments = load_alignments(["tests/data/titer_model/aa-seq_HA1.fasta"], ["HA1"])
    return alignments


@pytest.fixture
def titers():
    """Returns a path to a subset of titers.

    """
    return "tests/data/titer_model/h3n2_titers_subset.tsv"


class TestTiterTreeModel:
    def test_validate(self, tree, titers):
        model = TreeModel(tree, titers)
        model.prepare(training_fraction=0.5)
        model.train(method='nnl1reg')
        model_performance = model.validate()
        assert len(model_performance) > 0


class TestTiterSubstitutionModel:
    def test_validate(self, alignments, titers):
        model = SubstitutionModel(alignments, titers)
        model.prepare(training_fraction=0.5)
        model.train(method='nnl1reg')
        model_performance = model.validate()
        assert len(model_performance) > 0
