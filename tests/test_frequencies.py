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
    with open("tests/data/flu_seasonal_h3n2_ha_2y_tree.json") as fh:
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
