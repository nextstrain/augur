"""
Unit tests for frequency estimation
"""
import Bio
import json
import numpy as np
from pathlib import Path
import pytest
import sys
import os

# we assume (and assert) that this script is running from the tests/ directory
sys.path.append(str(Path(__file__).parent.parent.parent))

from augur.frequency_estimators import get_pivots, TreeKdeFrequencies, AlignmentKdeFrequencies
from augur.utils import json_to_tree

# Define regions to use for testing weighted frequencies.
REGIONS = [
    ('africa', 1.02),
    ('europe', 0.74),
    ('north_america', 0.54),
    ('china', 1.36),
    ('south_asia', 1.45),
    ('japan_korea', 0.2),
    ('oceania', 0.04),
    ('south_america', 0.41),
    ('southeast_asia', 0.62),
    ('west_asia', 0.75)
]

@pytest.fixture
def alignment():
    """Returns a multiple sequence alignment containing a small test set of H3N2
    sequences.
    """
    msa = Bio.AlignIO.read("tests/data/aa-seq_h3n2_ha_2y_HA1.fasta", "fasta")
    return msa

@pytest.fixture
def tree():
    """Returns an annotated Bio.Phylo tree.
    """
    with open("tests/data/flu_seasonal_h3n2_ha_3y_tree.json") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)
    return tree

#
# Test pivot calculations
#

def test_get_pivots_from_tree_only(tree):
    """Test pivot calculations.
    """
    # Define pivot frequency in months.
    pivot_frequency = 3
    observations = [tip.attr["num_date"] for tip in tree.get_terminals()]
    pivots = get_pivots(observations, pivot_frequency)
    assert isinstance(pivots, np.ndarray)

    # Floating point pivot values should be separated by the given number of
    # months divided by number of months in a year.
    assert pivots[1] - pivots[0] == pivot_frequency / 12.0

def test_get_pivots_from_start_and_end_date():
    """
    Test pivot calculation from a given start and end date instead of a given tree.
    """
    pivot_frequency = 3
    start_date = 2015.5
    end_date = 2018.5
    observations = []
    pivots = get_pivots(observations, pivot_frequency, start_date=start_date, end_date=end_date)
    assert isinstance(pivots, np.ndarray)
    assert pivots[1] - pivots[0] == pivot_frequency / 12.0
    assert pivots[0] == start_date
    assert pivots[-1] == end_date
    assert pivots[-1] >= end_date - pivot_frequency

#
# Test KDE frequency estimation for trees
#

class TestTreeKdeFrequencies(object):
    """Tests KDE-based frequency estimation methods for trees
    """
    def test_estimate(self, tree):
        """Test frequency estimation with default parameters.
        """
        kde_frequencies = TreeKdeFrequencies()
        frequencies = kde_frequencies.estimate(tree)
        assert hasattr(kde_frequencies, "pivots")
        assert np.around(kde_frequencies.pivots[1] - kde_frequencies.pivots[0], 2) == np.around(1 / 12.0, 2)
        assert hasattr(kde_frequencies, "frequencies")
        assert list(frequencies.values())[0].shape == kde_frequencies.pivots.shape

        # Frequencies should sum to 1 at all pivots.
        assert np.allclose(np.array(list(frequencies.values())).sum(axis=0), np.ones_like(kde_frequencies.pivots))

    def test_estimate_with_time_interval(self, tree):
        """Test frequency estimation with a given time interval.
        """
        start_date = 2015.5
        end_date = 2018.5
        kde_frequencies = TreeKdeFrequencies(
            start_date=start_date,
            end_date=end_date
        )
        frequencies = kde_frequencies.estimate(tree)
        assert hasattr(kde_frequencies, "pivots")
        assert kde_frequencies.pivots[0] == start_date
        assert hasattr(kde_frequencies, "frequencies")
        assert list(frequencies.values())[0].shape == kde_frequencies.pivots.shape

    def test_weighted_estimate(self, tree):
        """Test frequency estimation with weighted tips.
        """
        # Estimate weighted frequencies.
        weights = {region[0]: region[1] for region in REGIONS}
        kde_frequencies = TreeKdeFrequencies(
            weights=weights,
            weights_attribute="region"
        )
        frequencies = kde_frequencies.estimate(tree)
        assert hasattr(kde_frequencies, "pivots")
        assert hasattr(kde_frequencies, "frequencies")
        assert list(frequencies.values())[0].shape == kde_frequencies.pivots.shape

        # Frequencies should sum to 1 at all pivots.
        assert np.allclose(np.array(list(frequencies.values())).sum(axis=0), np.ones_like(kde_frequencies.pivots))

        # Estimate unweighted frequencies to compare with weighted frequencies.
        unweighted_kde_frequencies = TreeKdeFrequencies()
        unweighted_frequencies = unweighted_kde_frequencies.estimate(tree)

        # Any non-root node of the tree should have different frequencies with
        # or without weighting.
        clade_to_test = tree.root.clades[0]
        assert not np.array_equal(
            frequencies[clade_to_test.name],
            unweighted_frequencies[clade_to_test.name]
        )

    def test_weighted_estimate_with_unrepresented_weights(self, tree):
        """Test frequency estimation with weighted tips when any of the weight
        attributes is unrepresented.

        In this case, normalization of frequencies to the proportions
        represented by the weights should be followed by a second normalization
        to sum to 1.
        """
        # Drop all tips sampled from Africa from the tree. Despite dropping a
        # populous region, the estimated frequencies should still sum to 1
        # below.
        tips_from_africa = [
            tip
            for tip in tree.find_clades(terminal=True)
            if tip.attr["region"] == "africa"
        ]
        for tip in tips_from_africa:
            tree.prune(tip)

        # Estimate weighted frequencies.
        weights = {region[0]: region[1] for region in REGIONS}
        kde_frequencies = TreeKdeFrequencies(
            weights=weights,
            weights_attribute="region"
        )
        frequencies = kde_frequencies.estimate(tree)

        # Frequencies should sum to 1 at all pivots.
        assert np.allclose(np.array(list(frequencies.values())).sum(axis=0), np.ones_like(kde_frequencies.pivots))

    def test_only_tip_estimates(self, tree):
        """Test frequency estimation for only tips in a given tree.
        """
        # Estimate unweighted frequencies.
        kde_frequencies = TreeKdeFrequencies(
            include_internal_nodes=False
        )
        frequencies = kde_frequencies.estimate(tree)

        # Verify that all tips have frequency estimates and none of the internal nodes do.
        assert all([tip.name in frequencies
                    for tip in tree.get_terminals()])

        assert not any([node.clade in frequencies
                        for node in tree.get_nonterminals()])

        # Estimate weighted frequencies.
        weights = {region[0]: region[1] for region in REGIONS}
        kde_frequencies = TreeKdeFrequencies(
            weights=weights,
            weights_attribute="region",
            include_internal_nodes=False
        )
        frequencies = kde_frequencies.estimate(tree)

        # Verify that all tips have frequency estimates and none of the internal nodes do.
        assert all([tip.name in frequencies
                    for tip in tree.get_terminals()])

        assert not any([node.clade in frequencies
                        for node in tree.get_nonterminals()])

    def test_tip_and_internal_node_estimates(self, tree):
        """Test frequency estimation for tips and internal nodes in a given tree.
        """
        # Estimate unweighted frequencies.
        kde_frequencies = TreeKdeFrequencies(
            include_internal_nodes=True
        )
        frequencies = kde_frequencies.estimate(tree)

        # Verify that all tips and internal nodes have frequency estimates.
        assert all([tip.name in frequencies
                    for tip in tree.find_clades()])

    def test_censored_frequencies(self, tree):
        """Test estimation of frequencies where tips sampled beyond a given date are censored from the calculations.
        """
        max_date = 2017.0
        kde_frequencies = TreeKdeFrequencies(
            max_date=max_date
        )
        frequencies = kde_frequencies.estimate(tree)

        # Confirm that tips sampled after the max date have zero frequencies.
        assert all([frequencies[tip.name].sum() == 0
                    for tip in tree.get_terminals()
                    if tip.attr["num_date"] > max_date])

        # Confirm that one or more tips sampled before the max date have nonzero frequencies.
        assert any([frequencies[tip.name].sum() > 0
                    for tip in tree.get_terminals()
                    if tip.attr["num_date"] <= max_date])

    def test_node_filter(self, tree):
        """Test frequency estimation with specific nodes omitted by setting their
        frequencies to zero at all pivots.
        """
        # Filter nodes by region.
        regions = ["china"]
        kde_frequencies = TreeKdeFrequencies(
            node_filters={"region": regions}
        )
        frequencies = kde_frequencies.estimate(tree)

        # Verify that all tips have frequency estimates regardless of node
        # filter.
        assert all([tip.name in frequencies
                    for tip in tree.get_terminals()])

        # Verify that all tips from the requested region have non-zero frequencies.
        assert all([frequencies[tip.name].sum() > 0
                    for tip in tree.get_terminals()
                    if tip.attr["region"] in regions])

        # Verify that all tips not from the requested region have zero frequencies.
        assert all([frequencies[tip.name].sum() == 0
                    for tip in tree.get_terminals()
                    if tip.attr["region"] not in regions])

    def test_export_with_frequencies(self, tree):
        """Test frequencies export to JSON when frequencies have been estimated.
        """
        kde_frequencies = TreeKdeFrequencies()
        frequencies = kde_frequencies.estimate(tree)
        frequencies_json = kde_frequencies.to_json()

        assert "params" in frequencies_json
        assert kde_frequencies.pivot_frequency == frequencies_json["params"]["pivot_frequency"]
        assert kde_frequencies.start_date == frequencies_json["params"]["start_date"]
        assert kde_frequencies.end_date == frequencies_json["params"]["end_date"]
        assert "data" in frequencies_json
        assert "pivots" in frequencies_json["data"]
        assert "frequencies" in frequencies_json["data"]

    def test_export_without_frequencies(self):
        """Test frequencies export to JSON when frequencies have *not* been estimated.
        """
        kde_frequencies = TreeKdeFrequencies()
        frequencies_json = kde_frequencies.to_json()

        assert "params" in frequencies_json
        assert kde_frequencies.pivot_frequency == frequencies_json["params"]["pivot_frequency"]
        assert "node_filters" in frequencies_json["params"]
        assert "data" not in frequencies_json

    def test_import(self, tree, tmpdir):
        """Test import of frequencies JSON that was exported from a frequencies instance.
        """
        start_date = 2015.5
        end_date = 2018.5
        kde_frequencies = TreeKdeFrequencies(
            start_date=start_date,
            end_date=end_date
        )
        frequencies = kde_frequencies.estimate(tree)
        frequencies_json = kde_frequencies.to_json()

        # Try to dump exported JSON to disk.
        tmp_fh = tmpdir.mkdir("json").join("frequencies.json")
        fh = tmp_fh.open(mode="w")
        json.dump(frequencies_json, fh)
        fh.close()
        assert tmp_fh.check()

        # Import frequencies from existing tree and JSON.
        fh = tmp_fh.open()
        new_frequencies_json = json.load(fh)
        fh.close()
        new_kde_frequencies = TreeKdeFrequencies.from_json(new_frequencies_json)

        assert np.array_equal(
            kde_frequencies.pivots,
            new_kde_frequencies.pivots
        )

        # Get the first non-root key (root clade is number 0) and should be first in the sorted list of keys.
        key = sorted(kde_frequencies.frequencies.keys())[1]
        assert np.array_equal(
            kde_frequencies.frequencies[key],
            new_kde_frequencies.frequencies[key]
        )

    def test_import_without_frequencies(self):
        """Test import of frequencies JSON that was exported from a frequencies instance without frequency values.
        """
        kde_frequencies = TreeKdeFrequencies()
        frequencies_json = kde_frequencies.to_json()

        # Import frequencies from existing tree and JSON.
        new_kde_frequencies = TreeKdeFrequencies.from_json(frequencies_json)

        assert kde_frequencies.pivot_frequency == new_kde_frequencies.pivot_frequency
        assert not hasattr(new_kde_frequencies, "frequencies")

    def test_get_params(self, tree):
        """Test export of parameters used to create an instance.
        """
        initial_params = {
            "max_date": 2017.0,
            "start_date": 2015.5,
            "end_date": 2018.5
        }
        kde_frequencies = TreeKdeFrequencies(**initial_params)
        frequencies = kde_frequencies.estimate(tree)

        # Confirm that the exported parameters match the input.
        params = kde_frequencies.get_params()
        for param in initial_params:
            assert params[param] == initial_params[param]

#
# Test KDE frequency estimation for multiple sequence alignments
#

class TestAlignmentKdeFrequencies(object):
    """Tests KDE-based frequency estimation methods for multiple sequence alignments
    """
    def test_estimate(self, alignment):
        # Generate random dates for each sequence in the alignment.
        np.random.seed(1)
        observations = np.random.choice([2010.0, 2011.0], size=len(alignment)) + np.random.random(len(alignment))

        # Estimate frequencies for the alignment.
        kde_frequencies = AlignmentKdeFrequencies()
        frequencies = kde_frequencies.estimate(alignment, observations)

        assert hasattr(kde_frequencies, "pivots")
        assert np.around(kde_frequencies.pivots[1] - kde_frequencies.pivots[0], 2) == np.around(1 / 12.0, 2)
        assert hasattr(kde_frequencies, "frequencies")
        assert list(frequencies.values())[0].shape == kde_frequencies.pivots.shape

        # Find a position with frequencies estimated for multiple bases.
        position = list(frequencies.keys())[0].split(":")[0]
        position_frequencies = np.zeros_like(kde_frequencies.pivots)
        expected_position_frequencies = np.ones_like(kde_frequencies.pivots)
        for key in frequencies:
            if key.startswith("%s:" % position):
                position_frequencies += frequencies[key]

        # Confirm that the frequencies at this position sum to 1 at all timepoints.
        assert np.allclose(position_frequencies, expected_position_frequencies)
