"""
Tests for the `fitness_model` module.
"""
import Bio.Align.AlignInfo
import Bio.Phylo
import Bio.SeqIO
import datetime
import numpy as np
import pytest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from ..base.fitness_model import fitness_model

#
# Fixtures
#

# Set precalculated fitness model parameters which are the mean and standard
# deviation for the model.
MODEL_PARAMS = [1.0, 0.05]

@pytest.fixture
def simple_tree():
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
    sequences = (root, leaf_a, leaf_b)
    dates = (2012.5, 2013.25, 2014.8)
    index = 0
    for node in tree.find_clades(order="preorder"):
        node.aa = sequences[index]
        node.numdate = dates[index]
        index += 1

    return tree

@pytest.fixture
def real_tree(multiple_sequence_alignment):
    """Returns a tree built with FastTree from a small set of nucleotide sequences
    for H3N2.
    """
    # Load the tree.
    tree = Bio.Phylo.read("tests/fitness_model/H3N2_tree.newick", "newick")

    # Make a lookup table of name to sequence.
    sequences_by_name = dict([(alignment.name, str(alignment.seq))
                  for alignment in multiple_sequence_alignment])

    # Assign sequences to the tree.
    for node in tree.find_clades():
        if node.name is not None:
            node.sequence = np.fromstring(sequences_by_name[node.name], "S1")

            # Since sequence names look like "A/Singapore/TT0495/2017",
            # convert the last element to a floating point value for
            # simplicity.
            node.numdate = float(node.name.split("/")[-1])
        else:
            # Build a "dumb" consensus from the alignment for the
            # ancestral node and assign an arbitrary date in the
            # past.
            summary = Bio.Align.AlignInfo.SummaryInfo(multiple_sequence_alignment)
            node.sequence = np.fromstring(str(summary.dumb_consensus(threshold=0.5, ambiguous="N")), "S1")
            node.numdate = 2014.8

    return tree

@pytest.fixture
def simple_fitness_model(simple_tree):
    return fitness_model(
        tree=simple_tree,
        frequencies={},
        predictor_input=["random"],
        pivot_spacing=1.0 / 12,
        time_interval=(
            datetime.date(2015, 1, 1),
            datetime.date(2012, 1, 1)
        )
    )

@pytest.fixture
def real_fitness_model(real_tree, multiple_sequence_alignment):
    model = fitness_model(
        tree=real_tree,
        frequencies={},
        predictor_input=["random"],
        pivot_spacing=1.0 / 12,
        time_interval=(
            datetime.date(2017, 6, 1),
            datetime.date(2014, 6, 1)
        )
    )
    model.nuc_aln = multiple_sequence_alignment
    model.nuc_alphabet = 'ACGT-N'
    model.min_mutation_frequency = 0.01
    return model

@pytest.fixture
def precalculated_fitness_model(simple_tree):
    """Provides a simple fitness model with precalculated model parameters such that
    the model skips learning new parameters.
    """
    return fitness_model(
        tree=simple_tree,
        frequencies={},
        predictor_input={"random": MODEL_PARAMS},
        pivot_spacing=1.0 / 12,
        time_interval=(
            datetime.date(2015, 1, 1),
            datetime.date(2012, 1, 1)
        )
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
def multiple_sequence_alignment():
    """Returns a multiple sequence alignment containing a small test set of H3N2
    sequences.
    """
    msa = Bio.AlignIO.read("tests/fitness_model/H3N2_alignment.cleaned.fasta", "fasta")
    return msa

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

class TestFitnessModel(object):
    def test_prep_nodes(self, simple_fitness_model):
        assert not hasattr(simple_fitness_model, "nodes")
        assert not any([hasattr(node, "tips") for node in simple_fitness_model.tree.find_clades()])
        simple_fitness_model.prep_nodes()
        assert hasattr(simple_fitness_model, "nodes")
        assert hasattr(simple_fitness_model, "rootnode")
        assert hasattr(simple_fitness_model.rootnode, "pivots")
        assert all([hasattr(node, "tips") for node in simple_fitness_model.tree.find_clades()])

    def test_calc_node_frequencies(self, simple_fitness_model):
        simple_fitness_model.prep_nodes()
        assert not hasattr(simple_fitness_model, "freq_arrays")
        simple_fitness_model.calc_node_frequencies()
        assert hasattr(simple_fitness_model, "freq_arrays")
        assert len(simple_fitness_model.freq_arrays) > 0

    def test_calc_all_predictors(self, simple_fitness_model):
        simple_fitness_model.prep_nodes()
        simple_fitness_model.calc_node_frequencies()
        assert not hasattr(simple_fitness_model, "predictor_arrays")
        simple_fitness_model.calc_all_predictors()
        assert hasattr(simple_fitness_model, "predictor_arrays")
        assert len(simple_fitness_model.predictor_arrays) > 0

    def test_standardize_predictors(self, simple_fitness_model):
        simple_fitness_model.prep_nodes()
        simple_fitness_model.calc_node_frequencies()
        simple_fitness_model.calc_all_predictors()
        assert not hasattr(simple_fitness_model, "predictor_means")
        simple_fitness_model.standardize_predictors()
        assert hasattr(simple_fitness_model, "predictor_means")

    def test_select_clades_for_fitting(self, simple_fitness_model):
        simple_fitness_model.prep_nodes()
        simple_fitness_model.calc_node_frequencies()
        simple_fitness_model.calc_all_predictors()
        simple_fitness_model.standardize_predictors()
        assert not hasattr(simple_fitness_model, "fit_clades")
        simple_fitness_model.select_clades_for_fitting()
        assert hasattr(simple_fitness_model, "fit_clades")
        assert len(simple_fitness_model.fit_clades) > 0

    def test_learn_parameters(self, real_fitness_model):
        real_fitness_model.prep_nodes()
        real_fitness_model.calc_node_frequencies()
        real_fitness_model.calc_all_predictors()
        real_fitness_model.standardize_predictors()
        real_fitness_model.select_clades_for_fitting()
        assert not hasattr(real_fitness_model, "last_fit")
        real_fitness_model.learn_parameters(niter=1, fit_func="clade")
        assert hasattr(real_fitness_model, "last_fit")

    def test_assign_fitness(self, real_fitness_model):
        real_fitness_model.prep_nodes()
        real_fitness_model.calc_node_frequencies()
        real_fitness_model.calc_all_predictors()
        real_fitness_model.standardize_predictors()
        real_fitness_model.select_clades_for_fitting()
        real_fitness_model.learn_parameters(niter=1, fit_func="clade")
        assert not any([hasattr(node, "fitness") for node in real_fitness_model.tree.find_clades()])
        real_fitness_model.assign_fitness()
        assert all([hasattr(node, "fitness") for node in real_fitness_model.tree.find_clades()])

    def test_assign_fitness_with_precalculated_params(self, precalculated_fitness_model):
        # The fitness model should have model parameters assigned by the user.
        assert np.array_equal(precalculated_fitness_model.model_params, np.array([MODEL_PARAMS[0]]))
        precalculated_fitness_model.predict()

        # After prediction, the model parameters should be unchanged as the
        # learning step should be skipped.
        assert np.array_equal(precalculated_fitness_model.model_params, np.array([MODEL_PARAMS[0]]))

        # Recalculate fitness model parameters which should be different from those given.
        precalculated_fitness_model.learn_parameters(niter=1, fit_func="clade")
        assert not np.array_equal(precalculated_fitness_model.model_params, np.array([MODEL_PARAMS[0]]))
