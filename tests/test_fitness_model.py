"""
Tests for the `fitness_model` module.
"""
import Bio.Align.AlignInfo
import Bio.Phylo
import Bio.SeqIO
import dendropy
import pytest

try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO

#
# Fixtures
#

@pytest.fixture
def simple_tree():
	"""Returns a tree with three sequences: a root and two direct descendents with
	one modification each.
	"""
	# Build simple tree.
	tree = dendropy.Tree(stream=StringIO("(A,B);"), schema="newick")

	# Build sequences for tree nodes. One leaf has a Koel and epitope site
	# mutation. The other leaf has a signal peptide mutation.
	root = sequence()
	leaf_a = modify_sequence_at_site(root, 145 + 16 - 1)
	leaf_b = modify_sequence_at_site(root, 14)

	# Assign sequences to nodes.
	sequences = (root, leaf_a, leaf_b)
	dates = (2012.5, 2013.25, 2014.8)
	index = 0
	for node in tree.preorder_node_iter():
		node.aa = sequences[index]
		node.num_date = dates[index]
		index += 1

	return tree

@pytest.fixture
def real_tree(multiple_sequence_alignment):
	"""Returns a tree built with FastTree from a small set of nucleotide sequences
	for H3N2.
	"""
	# Load the tree.
	tree = dendropy.Tree.get_from_path("tests/data/H3N2_tree.newick", "newick")

	# Make a lookup table of name to sequence.
	sequences_by_name = dict([(alignment.name, str(alignment.seq))
				  for alignment in multiple_sequence_alignment])

	# Assign sequences to the tree.
	for node in tree.nodes():
		if node.taxon is not None:
			node.seq = sequences_by_name[node.taxon.label]

			# Since sequence names look like "A/Singapore/TT0495/2017",
			# convert the last element to a floating point value for
			# simplicity.
			node.num_date = float(node.taxon.label.split("/")[-1])
		else:
			# Build a "dumb" consensus from the alignment for the
			# ancestral node and assign an arbitrary date in the
			# past.
			summary = Bio.Align.AlignInfo.SummaryInfo(multiple_sequence_alignment)
			node.seq = str(summary.dumb_consensus(threshold=0.5, ambiguous="N"))
			node.num_date = 2014.8

	return tree

@pytest.fixture
def fitness_model(simple_tree):
	from src.fitness_model import fitness_model
	return fitness_model(
		predictor_input=["ep"],
		pivots_per_year=12,
		time_interval=(2012.0, 2015.0),
		tree=simple_tree
	)

@pytest.fixture
def real_fitness_model(real_tree, multiple_sequence_alignment):
	from src.fitness_model import fitness_model
	model = fitness_model(
		predictor_input=["ep"],
		pivots_per_year=12,
		time_interval=(2014.5, 2017.5),
		tree=real_tree
	)
	model.nuc_aln = multiple_sequence_alignment
	model.nuc_alphabet = 'ACGT-N'
	model.min_mutation_frequency = 0.01
	return model

@pytest.fixture
def sequence():
	"""Returns an amino acid sequence for an ancestral H3N2 virus (Hong Kong 1968).
	"""
	with open("tests/data/AAK51718.fasta", "r") as handle:
		record = list(Bio.SeqIO.parse(handle, "fasta"))[0]

	aa = str(record.seq)
	return aa

@pytest.fixture
def multiple_sequence_alignment():
	"""Returns a multiple sequence alignment containing a small test set of H3N2
	sequences.
	"""
	msa = Bio.AlignIO.read("tests/data/H3N2_alignment.cleaned.fasta", "fasta")
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
	def test_prep_nodes(self, fitness_model):
		assert not hasattr(fitness_model, "nodes")
		fitness_model.prep_nodes()
		assert hasattr(fitness_model, "nodes")
		assert hasattr(fitness_model, "rootnode")
		assert hasattr(fitness_model.rootnode, "pivots")

	def test_calc_node_frequencies(self, fitness_model):
		fitness_model.prep_nodes()
		assert not hasattr(fitness_model, "freq_arrays")
		fitness_model.calc_node_frequencies()
		assert hasattr(fitness_model, "freq_arrays")
		assert len(fitness_model.freq_arrays) > 0

	def test_calc_all_predictors(self, fitness_model):
		fitness_model.prep_nodes()
		fitness_model.calc_node_frequencies()
		assert not hasattr(fitness_model, "predictor_arrays")
		fitness_model.calc_all_predictors()
		assert hasattr(fitness_model, "predictor_arrays")
		assert len(fitness_model.predictor_arrays) > 0

	def test_standardize_predictors(self, fitness_model):
		fitness_model.prep_nodes()
		fitness_model.calc_node_frequencies()
		fitness_model.calc_all_predictors()
		assert not hasattr(fitness_model, "predictor_means")
		fitness_model.standardize_predictors()
		assert hasattr(fitness_model, "predictor_means")

	def test_select_clades_for_fitting(self, fitness_model):
		fitness_model.prep_nodes()
		fitness_model.calc_node_frequencies()
		fitness_model.calc_all_predictors()
		fitness_model.standardize_predictors()
		assert not hasattr(fitness_model, "fit_clades")
		fitness_model.select_clades_for_fitting()
		assert hasattr(fitness_model, "fit_clades")
		assert len(fitness_model.fit_clades) > 0

	def test_prep_af(self, real_fitness_model):
		real_fitness_model.prep_nodes()
		real_fitness_model.calc_node_frequencies()
		real_fitness_model.calc_all_predictors()
		real_fitness_model.standardize_predictors()
		real_fitness_model.select_clades_for_fitting()
		assert not hasattr(real_fitness_model, "af")
		real_fitness_model.prep_af()
		assert hasattr(real_fitness_model, "af")
		assert len(real_fitness_model.af) > 0

	def test_learn_parameters(self, real_fitness_model):
		real_fitness_model.prep_nodes()
		real_fitness_model.calc_node_frequencies()
		real_fitness_model.calc_all_predictors()
		real_fitness_model.standardize_predictors()
		real_fitness_model.select_clades_for_fitting()
		real_fitness_model.prep_af()
		assert not hasattr(real_fitness_model, "last_fit")
		real_fitness_model.learn_parameters(niter=1, fit_func="clade")
		assert hasattr(real_fitness_model, "last_fit")

	def test_assign_fitness(self, real_fitness_model):
		real_fitness_model.prep_nodes()
		real_fitness_model.calc_node_frequencies()
		real_fitness_model.calc_all_predictors()
		real_fitness_model.standardize_predictors()
		real_fitness_model.select_clades_for_fitting()
		real_fitness_model.prep_af()
		real_fitness_model.learn_parameters(niter=1, fit_func="clade")
		assert not any([hasattr(node, "fitness") for node in real_fitness_model.tree.nodes()])
		real_fitness_model.assign_fitness()
		assert all([hasattr(node, "fitness") for node in real_fitness_model.tree.nodes()])
