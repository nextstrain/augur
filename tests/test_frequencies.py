"""
Unit tests for frequency estimation
"""
import json
import numpy as np
import pytest

from ..base.frequencies import KdeFrequencies
from ..base.io_util import json_to_tree

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
def tree():
    """Returns an annotated Bio.Phylo tree.
    """
    with open("tests/data/flu_seasonal_h3n2_ha_3y_tree.json") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)
    return tree


class TestKdeFrequencies(object):
    """Tests KDE-based frequency estimation methods
    """
    def test_calculate_pivots(self, tree):
        """Test pivot calculations.
        """
        pivot_frequency = 0.25
        pivots = KdeFrequencies.calculate_pivots(tree, pivot_frequency)
        assert isinstance(pivots, np.ndarray)
        assert pivots[1] - pivots[0] == pivot_frequency

    def test_estimate(self, tree):
        """Test frequency estimation with default parameters.
        """
        kde_frequencies = KdeFrequencies()
        frequencies = kde_frequencies.estimate(tree)
        assert "global" in frequencies
        assert hasattr(kde_frequencies, "pivots")
        assert hasattr(kde_frequencies, "frequencies")
        assert frequencies["global"].values()[0].shape == kde_frequencies.pivots.shape

    def test_weighted_estimate(self, tree):
        """Test frequency estimation with weighted tips.
        """
        # Estimate weighted frequencies.
        weights = {region[0]: region[1] for region in REGIONS}
        kde_frequencies = KdeFrequencies(
            weights=weights,
            weights_attribute="region"
        )
        frequencies = kde_frequencies.estimate(tree)
        assert "global" in frequencies
        assert hasattr(kde_frequencies, "pivots")
        assert hasattr(kde_frequencies, "frequencies")
        assert frequencies["global"].values()[0].shape == kde_frequencies.pivots.shape

        # Estimate unweighted frequencies to compare with weighted frequencies.
        unweighted_kde_frequencies = KdeFrequencies()
        unweighted_frequencies = unweighted_kde_frequencies.estimate(tree)

        # The any non-root node of the tree should have different frequencies with or without weighting.
        assert not np.array_equal(
            frequencies["global"][1],
            unweighted_frequencies["global"][1]
        )

    def test_only_tip_estimates(self, tree):
        """Test frequency estimation for only tips in a given tree.
        """
        kde_frequencies = KdeFrequencies(
            include_internal_nodes=False
        )
        frequencies = kde_frequencies.estimate(tree)

        # Verify that all tips have frequency estimates and none of the internal nodes do.
        assert all([tip.clade in frequencies["global"]
                    for tip in tree.get_terminals()])

        assert not any([node.clade in frequencies["global"]
                        for node in tree.get_nonterminals()])

    def test_censored_frequencies(self, tree):
        """Test estimation of frequencies where tips sampled beyond a given date are censored from the calculations.
        """
        max_date = 2017.0
        kde_frequencies = KdeFrequencies(
            max_date=max_date
        )
        frequencies = kde_frequencies.estimate(tree)

        # Confirm that tips sampled after the max date have zero frequencies.
        assert all([frequencies["global"][tip.clade].sum() == 0
                    for tip in tree.get_terminals()
                    if tip.attr["num_date"] > max_date])

        # Confirm that one or more tips sampled before the max date have nonzero frequencies.
        assert any([frequencies["global"][tip.clade].sum() > 0
                    for tip in tree.get_terminals()
                    if tip.attr["num_date"] <= max_date])

    def test_export_with_frequencies(self, tree, tmpdir):
        """Test frequencies export to JSON when frequencies have been estimated.
        """
        kde_frequencies = KdeFrequencies()
        frequencies = kde_frequencies.estimate(tree)
        frequencies_json = kde_frequencies.to_json()

        assert "params" in frequencies_json
        assert kde_frequencies.pivot_frequency == frequencies_json["params"]["pivot_frequency"]
        assert "data" in frequencies_json
        assert "pivots" in frequencies_json["data"]
        assert "frequencies" in frequencies_json["data"]
        assert "global" in frequencies_json["data"]["frequencies"]

        # Try to dump exported JSON to disk.
        fh = tmpdir.mkdir("json").join("frequencies.json")
        json.dump(frequencies_json, fh)

    def test_export_without_frequencies(self):
        """Test frequencies export to JSON when frequencies have *not* been estimated.
        """
        kde_frequencies = KdeFrequencies()
        frequencies_json = kde_frequencies.to_json()

        assert "params" in frequencies_json
        assert kde_frequencies.pivot_frequency == frequencies_json["params"]["pivot_frequency"]
        assert "data" not in frequencies_json
