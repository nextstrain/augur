import pytest
import pandas as pd

@pytest.fixture
def geolocation_rules():
    """
    Returns the geolocation rules as a Pandas DataFrame
    """
    geolocation_rules_file = "augur/data/geolocation_rules.tsv"
    return pd.read_csv(
        geolocation_rules_file,
        sep="\t",
        names=['raw', 'annotations'],
        comment='#')


def test_geolocation_rules_are_not_duplicated(geolocation_rules):
    """
    Check for duplicated raw values.
    We want to avoid overwriting rules within our own rules.
    """
    duplicate_rule = geolocation_rules.loc[geolocation_rules.raw.duplicated()]
    assert len(duplicate_rule) == 0


def test_geolocation_rules_are_changing_values(geolocation_rules):
    """
    Check for rules where the raw values and annotations are the same.
    These rules are unnecessary as the values will be unchanged.
    """
    unchanged_rule = geolocation_rules.loc[geolocation_rules.raw == geolocation_rules.annotations]
    assert len(unchanged_rule) == 0


def test_geolocation_rules_are_not_cyclic(geolocation_rules):
    """
    Check for rules where the annotations exist in the raw values and
    the raw values exist in the annotations since this is the simplest
    way to flag cyclic rules.

    It won't catch all cyclic rules because wildcards might mask them,
    but it's better than nothing!
    """
    # Filter rules to ones that are changing values to make sure the same rules
    # don't get flagged by this test and `test_geolocation_rules_are_changing_values`
    changed_rules = geolocation_rules.loc[geolocation_rules.raw != geolocation_rules.annotations]
    cyclic_rules = changed_rules.loc[
        changed_rules.annotations.isin(changed_rules.raw) &
        changed_rules.raw.isin(changed_rules.annotations)
    ]
    assert len(cyclic_rules) == 0
