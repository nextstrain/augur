"""
Unit tests for frequency estimation
"""
import Bio
import json
import numpy as np
from pathlib import Path
import pytest
import sys

# we assume (and assert) that this script is running from the tests/ directory
sys.path.append(str(Path(__file__).parent.parent.parent))

from augur.dates import numeric_date
from augur.frequency_estimators import float_to_datestring, get_pivots, TreeKdeFrequencies, AlignmentKdeFrequencies, TreeKdeFrequenciesError
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

    # Floating point pivot values should be roughly separated by the given
    # number of months divided by number of months in a year. Numeric dates from
    # pivots calculate decimal fractions as the proportion of days in the year,
    # instead of months in the year (e.g., the fraction for the first three
    # months of a non-leap year is (31 + 28 + 31) / 365 or 0.246575 instead of 3
    # / 12 or 0.25). As a result of this difference in how we calculate
    # fractions, we need to round to pivots to 2 decimal places which is the
    # precision of month-based fractions.
    assert np.round(pivots[1] - pivots[0], 2) == (pivot_frequency / 12.0)

def test_get_pivots_from_start_and_end_date():
    """
    Test pivot calculation from a given start and end date instead of a given tree
    First pivot is the first day of the month immediately following start_date
    Last pivot is the first day of the month immediately preceding end_date
    Current logic converts numeric date 2015.5 to 2015-07-02, hence using 2015.49
    """
    pivot_frequency = 1
    start_date = 2015.49
    end_date = 2018.5
    observations = []
    pivots = get_pivots(observations, pivot_frequency, start_date=start_date, end_date=end_date)
    assert isinstance(pivots, np.ndarray)
    assert np.round( 12 * (pivots[1] - pivots[0]) ) == pivot_frequency
    assert pivots[0] == 2015.5
    assert pivots[-1] == 2018.5
    assert pivots[-1] >= end_date - pivot_frequency

def test_get_pivots_by_months():
    """Get pivots where intervals are defined by months.
    """
    pivots = get_pivots(observations=[], pivot_interval=1, start_date=2015.0, end_date=2016.0, pivot_interval_units="months")
    # Pivots should include all 12 months of the year plus the month represented
    # by the end date because pivots represent time slices through data and we
    # want to evaluate frequencies at the start and end time slices.
    assert len(pivots) == 13
    assert float_to_datestring(pivots[-1]) == "2016-01-01"

def test_get_pivots_by_months_with_realistic_start_end_dates():
    """Get pivots where intervals are defined by months and realistic start and end dates.
    """
    # 6 years (72 months) of pivots with 3-month intervals should produce 24 + 1
    # pivots (4 pivots per year plus the pivot for the beginning of the year
    # associated with the end date).
    pivots = get_pivots(
        observations=[],
        pivot_interval=3,
        start_date=numeric_date("2017-01-06"),
        end_date=numeric_date("2023-01-06"),
        pivot_interval_units="months"
    )
    assert len(pivots) == 25
    assert float_to_datestring(pivots[-1]) == "2023-01-06"

def test_get_pivots_by_weeks():
    """Get pivots where intervals are defined as weeks instead of months.
    """
    # As with monthly pivots, weekly pivots should include the first and last
    # values in the range. So, a 1-year interval will represent 52 weekly pivots
    # from the beginning plus the first pivot from the next year.
    pivots = get_pivots(observations=[], pivot_interval=1, start_date=2015.0, end_date=2016.0, pivot_interval_units="weeks")
    assert len(pivots) == 53

    pivots = get_pivots(
        observations=[],
        pivot_interval=1,
        start_date=numeric_date("2022-01-06"),
        end_date=numeric_date("2023-01-06"),
        pivot_interval_units="weeks"
    )
    assert len(pivots) == 53


def test_get_pivots_by_invalid_unit():
    with pytest.raises(ValueError, match=r".*invalid_unit.*is not supported.*"):
        pivots = get_pivots(observations=[], pivot_interval=1, start_date=2015.0, end_date=2016.0, pivot_interval_units="invalid_unit")


@pytest.mark.parametrize(
    "start, end, expected_pivots",
    [
        (
            "2022-01-01",
            "2022-04-01",
            ("2022-01-01", "2022-02-01", "2022-03-01", "2022-04-01")
        ),
        (
            "2022-01-31",
            "2022-03-31",
            ("2022-01-31", "2022-02-28", "2022-03-31")
        ),
        # Note that Jan 31 to Apr 30 gives the same amount of pivot points as
        # Jan 31 to Mar 31.
        (
            "2022-01-31",
            "2022-04-30",
            ("2022-02-28", "2022-03-30", "2022-04-30")
        ),
        # However, in practice, the interval is more likely to be Jan 30 to Apr
        # 30 as long as the start date is calculated relative to the end date
        # (i.e. start date = 3 months before Apr 30 = Jan 30).
        # That interval includes an additional pivot point as expected.
        (
            "2022-01-30",
            "2022-04-30",
            ("2022-01-30", "2022-02-28", "2022-03-30", "2022-04-30")
        ),
    ]
)
def test_get_pivots_on_month_boundaries(start, end, expected_pivots):
    """Get pivots where the start/end dates are on month boundaries.
    """
    pivots = get_pivots(
        observations=[],
        pivot_interval=1,
        start_date=numeric_date(start),
        end_date=numeric_date(end),
        pivot_interval_units="months"
    )
    assert len(pivots) == len(expected_pivots)
    assert np.allclose(pivots, [numeric_date(date) for date in expected_pivots], rtol=0, atol=1e-4)

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
        assert kde_frequencies.pivots[0] == 2015.5
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

        # Estimate weighted frequencies such that all weighted attributes are
        # missing. This should raise an exception because none of the tips will
        # match any of the weights and the weighting of frequencies will be
        # impossible.
        weights = {"fake_region_1": 1.0, "fake_region_2": 2.0}
        kde_frequencies = TreeKdeFrequencies(
            weights=weights,
            weights_attribute="region"
        )

        with pytest.raises(TreeKdeFrequenciesError):
            frequencies = kde_frequencies.estimate(tree)

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
        pivots = kde_frequencies.pivots

        # Pivot decimal fractions represent fractions of a year by days, so
        # comparisons with fractions of a year by months need to take into
        # account rounding error between these calculations. In this test, we
        # allow pivots to be spaced apart by a value that is "close" to 1 / 12
        # months by the rounding error of the month-based precision.
        assert np.isclose((pivots[-1] - pivots[-2]), 1 / 12.0, atol=0.005)
        assert hasattr(kde_frequencies, "frequencies")
        assert list(frequencies.values())[0].shape == pivots.shape

        # Find a position with frequencies estimated for multiple bases.
        position = list(frequencies.keys())[0].split(":")[0]
        position_frequencies = np.zeros_like(kde_frequencies.pivots)
        expected_position_frequencies = np.ones_like(kde_frequencies.pivots)
        for key in frequencies:
            if key.startswith("%s:" % position):
                position_frequencies += frequencies[key]

        # Confirm that the frequencies at this position sum to 1 at all timepoints.
        assert np.allclose(position_frequencies, expected_position_frequencies)
