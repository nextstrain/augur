from collections import defaultdict
from io import StringIO

import Bio.Phylo
import pytest

from augur.distance import get_distances_to_root
from augur.errors import AugurError


@pytest.fixture
def tree():
    return Bio.Phylo.read(StringIO("(A:1,(B:1,C:1)internal:1)root;"), "newick")


@pytest.fixture
def distance_map():
    return {"default": 1, "map": {}}


class TestGetDistancesToRoot:
    def test_with_root_sequence(self, tree, distance_map):
        sequences_by_node_and_gene = {
            "root": {"nuc": "ACGT"},
            "internal": {"nuc": "ACGA"},
            "A": {"nuc": "ACGT"},
            "B": {"nuc": "ACTA"},
            "C": {"nuc": "TCGA"},
        }

        distances = get_distances_to_root(tree, sequences_by_node_and_gene, distance_map)

        assert distances == {
            "root": 0.0,
            "internal": 1.0,
            "A": 0.0,
            "B": 2.0,
            "C": 2.0,
        }

    def test_without_root_sequence(self, tree, distance_map):
        # Only tips are present, as is the case when the alignment does not
        # include ancestral sequences for internal nodes.
        sequences_by_node_and_gene = {
            "A": {"nuc": "ACGT"},
            "B": {"nuc": "ACTA"},
            "C": {"nuc": "TCGA"},
        }

        with pytest.raises(AugurError, match="Could not find a sequence for the root node 'root'"):
            get_distances_to_root(tree, sequences_by_node_and_gene, distance_map)

    def test_without_root_sequence_in_defaultdict(self, tree, distance_map):
        # `augur distance` indexes sequences in a defaultdict, which previously
        # masked the missing root as an empty sequence and reported all
        # distances as zero.
        sequences_by_node_and_gene = defaultdict(dict)
        sequences_by_node_and_gene["A"]["nuc"] = "ACGT"
        sequences_by_node_and_gene["B"]["nuc"] = "ACTA"

        with pytest.raises(AugurError, match="Could not find a sequence for the root node 'root'"):
            get_distances_to_root(tree, sequences_by_node_and_gene, distance_map)
