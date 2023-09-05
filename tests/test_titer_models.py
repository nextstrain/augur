
from augur.titer_model import TiterCollection


def test_titer_collection():
    # Confirm that titers load from a file path.
    titers = TiterCollection("tests/data/titer_model/h3n2_titers_subset.tsv")

    # Confirm that all distinct test and reference strains have been counted.
    assert len(titers.strains) == 13
