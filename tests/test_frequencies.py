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
    def test_calculate_pivots_from_tree_only(self, tree):
        """Test pivot calculations.
        """
        pivot_frequency = 0.25
        pivots = KdeFrequencies.calculate_pivots(pivot_frequency, tree=tree)
        assert isinstance(pivots, np.ndarray)
        assert pivots[1] - pivots[0] == pivot_frequency

    def test_calculate_pivots_from_start_and_end_date(self):
        """
        Test pivot calculation from a given start and end date instead of a given tree.
        """
        pivot_frequency = 0.25
        start_date = 2015.5
        end_date = 2018.5
        pivots = KdeFrequencies.calculate_pivots(pivot_frequency, start_date=start_date, end_date=end_date)
        assert isinstance(pivots, np.ndarray)
        assert pivots[1] - pivots[0] == pivot_frequency
        assert pivots[0] == start_date
        assert pivots[-1] != end_date
        assert pivots[-1] >= end_date - pivot_frequency

    def test_estimate(self, tree):
        """Test frequency estimation with default parameters.
        """
        kde_frequencies = KdeFrequencies()
        frequencies = kde_frequencies.estimate(tree)
        assert "global" in frequencies
        assert hasattr(kde_frequencies, "pivots")
        assert hasattr(kde_frequencies, "frequencies")
        assert frequencies["global"].values()[0].shape == kde_frequencies.pivots.shape

    def test_estimate_with_time_interval(self, tree):
        """Test frequency estimation with a given time interval.
        """
        start_date = 2015.5
        end_date = 2018.5
        kde_frequencies = KdeFrequencies(
            start_date=start_date,
            end_date=end_date
        )
        frequencies = kde_frequencies.estimate(tree)
        assert "global" in frequencies
        assert hasattr(kde_frequencies, "pivots")
        assert kde_frequencies.pivots[0] == start_date
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
        assert kde_frequencies.start_date == frequencies_json["params"]["start_date"]
        assert kde_frequencies.end_date == frequencies_json["params"]["end_date"]
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

    def test_import(self, tree, tmpdir):
        """Test import of frequencies JSON that was exported from a frequencies instance.
        """
        start_date = 2015.5
        end_date = 2018.5
        kde_frequencies = KdeFrequencies(
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
        new_kde_frequencies = KdeFrequencies.from_json(new_frequencies_json)

        assert np.array_equal(
            kde_frequencies.pivots,
            new_kde_frequencies.pivots
        )

        # Get the first non-root key (root clade is number 0) and should be first in the sorted list of keys.
        key = sorted(kde_frequencies.frequencies["global"].keys())[1]
        assert np.array_equal(
            kde_frequencies.frequencies["global"][key],
            new_kde_frequencies.frequencies["global"][key]
        )

    def test_import_without_frequencies(self):
        """Test import of frequencies JSON that was exported from a frequencies instance without frequency values.
        """
        kde_frequencies = KdeFrequencies()
        frequencies_json = kde_frequencies.to_json()

        # Import frequencies from existing tree and JSON.
        new_kde_frequencies = KdeFrequencies.from_json(frequencies_json)

        assert kde_frequencies.pivot_frequency == new_kde_frequencies.pivot_frequency
        assert not hasattr(new_kde_frequencies, "frequencies")
