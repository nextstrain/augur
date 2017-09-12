"""
Tests for the `fitness_predictors` module.
"""
import Bio.Phylo
import Bio.SeqIO
import pytest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

#
# Fixtures
#

@pytest.fixture
def fitness_predictor():
    from ..base.fitness_predictors import fitness_predictors
    return fitness_predictors(
        epitope_mask_version="wolf",
        tolerance_mask_version="ha1"
    )

@pytest.fixture
def sequence():
    """Returns an amino acid sequence for an ancestral H3N2 virus (Hong Kong 1968).
    """
    with open("tests/fitness_model/AAK51718.fasta", "r") as handle:
        record = list(Bio.SeqIO.parse(handle, "fasta"))[0]

    aa = str(record.seq)
    return aa

@pytest.fixture
def tree():
    """Returns a tree with three sequences: a root and two direct descendents with
    one modification each.
    """
    # Build simple tree.
    tree = Bio.Phylo.read(StringIO("(A,B);"), "newick")

    # Build sequences for tree nodes. One leaf has a Koel and epitope site
    # mutation. The other leaf has a signal peptide mutation.
    root = sequence()
    leaf_a = modify_sequence_at_site(root, 145 + 16 - 1)
    leaf_b = modify_sequence_at_site(root, 14)

    # Assign sequences to nodes.
    sequences = [root, leaf_a, leaf_b]
    index = 0
    for node in tree.find_clades(order="preorder"):
        node.aa = sequences[index]
        index += 1

    return tree

#
# Utility functions
#

def modify_sequence_at_site(sequence, site):
    """Returns the given sequence with a modified base at the given site.
    """
    other_sequence_list = list(sequence)
    other_sequence_list[site] = "Z"
    return "".join(other_sequence_list)

#
# Tests
#

class TestFitnessPredictors(object):
    def test_epitope_sites(self, fitness_predictor, sequence):
        eps = fitness_predictor.epitope_sites(sequence)
        assert len(eps) == len(fitness_predictor.epitope_mask.replace("0", ""))

    def test_nonepitope_sites(self, fitness_predictor, sequence):
        non_epitope_sites = fitness_predictor.nonepitope_sites(sequence)

        # Non-epitope sites should be zero sites in the tolerance mask.
        assert len(non_epitope_sites) == len(fitness_predictor.tolerance_mask.replace("1", ""))

    def test_receptor_binding_sites(self, fitness_predictor, sequence):
        rbs = fitness_predictor.receptor_binding_sites(sequence)

        # There should be 7 Koel recepter binding sites.
        assert len(rbs) == 7

        # The Hong Kong 1968 ancestral sequence at the second Koel site was a
        # "T". See Koel et al. 2013 for more.
        assert rbs[1] == "T"

    def test_epitope_distance(self, fitness_predictor, sequence):
        # Create a copy of the given sequence.
        other_sequence = sequence
        assert fitness_predictor.epitope_distance(sequence, other_sequence) == 0

        # Modify the given sequence with a change at an epitope
        # site. Replacement is a non-standard amino acid for testing purposes.
        epitope_site_index = fitness_predictor.epitope_mask.index("1")
        other_sequence = modify_sequence_at_site(other_sequence, epitope_site_index)

        assert sequence != other_sequence
        assert fitness_predictor.epitope_distance(sequence, other_sequence) == 1

    def test_fast_epitope_distance(self, fitness_predictor, sequence):
        # Create a copy of the given sequence.
        other_sequence = sequence
        assert fitness_predictor.fast_epitope_distance(sequence, other_sequence) == 0

        # Modify the given sequence with a change at an epitope
        # site. Replacement is a non-standard amino acid for testing purposes.
        epitope_site_index = fitness_predictor.epitope_mask.index("1")
        other_sequence = modify_sequence_at_site(other_sequence, epitope_site_index)

        assert sequence != other_sequence
        assert fitness_predictor.fast_epitope_distance(sequence, other_sequence) == 1
        assert fitness_predictor.fast_epitope_distance(sequence, other_sequence) == fitness_predictor.epitope_distance(sequence, other_sequence)

    def test_nonepitope_distance(self, fitness_predictor, sequence):
        # Create a copy of the given sequence.
        other_sequence = sequence
        assert fitness_predictor.nonepitope_distance(sequence, other_sequence) == 0

        # Modify the given sequence with a change at a non-epitope
        # site. Replacement is a non-standard amino acid for testing purposes.
        nonepitope_site_index = fitness_predictor.tolerance_mask.index("0")
        other_sequence = modify_sequence_at_site(other_sequence, nonepitope_site_index)

        assert sequence != other_sequence
        assert fitness_predictor.nonepitope_distance(sequence, other_sequence) == 1

    def test_rbs_distance(self, fitness_predictor, sequence):
        # Create a copy of the given sequence.
        other_sequence = sequence
        assert fitness_predictor.rbs_distance(sequence, other_sequence) == 0

        # Modify the given sequence with a change at a Koel et al. 2013 site
        # (HA1 coordinates plus signal peptide length minus 1 for python
        # indexing). Replacement is a non-standard amino acid for testing
        # purposes.
        koel_site_index = 145 + 16 - 1
        other_sequence = modify_sequence_at_site(other_sequence, koel_site_index)

        assert sequence != other_sequence
        assert fitness_predictor.rbs_distance(sequence, other_sequence) == 1

    def test_calc_epitope_distance(self, fitness_predictor, tree):
        attr = "ep"
        nodes = [node for node in tree.find_clades(order="preorder")]
        for node in nodes:
            assert hasattr(node, "aa")
            assert not hasattr(node, attr)

        fitness_predictor.calc_epitope_distance(tree, attr=attr)

        # All nodes should have an epitope distance from the root.
        for node in nodes:
            assert hasattr(node, attr)

        # Only one node should have a non-zero epitope distance.
        assert getattr(nodes[0], attr) == 0
        assert getattr(nodes[1], attr) == 1
        assert getattr(nodes[2], attr) == 0

    def test_calc_null_predictor(self, fitness_predictor, tree):
        attr = "null"
        nodes = [node for node in tree.find_clades(order="preorder")]
        for node in nodes:
            assert not hasattr(node, attr)

        fitness_predictor.calc_null_predictor(tree, attr=attr)

        # All nodes should have a zero-valued predictor.
        for node in nodes:
            assert getattr(node, attr) == 0

    def test_calc_random_predictor(self, fitness_predictor, tree):
        attr = "random"
        nodes = [node for node in tree.find_clades(order="preorder")]
        for node in nodes:
            assert not hasattr(node, attr)

        fitness_predictor.calc_random_predictor(tree, attr=attr)

        # All nodes should have a random value between 0 and 1.
        for node in nodes:
            assert getattr(node, attr) >= 0 and getattr(node, attr) < 1
