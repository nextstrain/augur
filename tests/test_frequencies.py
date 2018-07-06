"""
Unit tests for frequency estimation
"""
import json
import numpy as np
import pytest

from ..base.frequencies import KdeFrequencies
from ..base.io_util import json_to_tree


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
        assert frequencies["global"].values()[0].shape == kde_frequencies.pivots.shape
