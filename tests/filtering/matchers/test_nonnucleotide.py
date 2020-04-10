import argparse

from augur.filtering.matchers.nonnucleotide import Nonnucleotide

import pytest


@pytest.fixture
def all_sequences(sequence_factory):
    return {
        sequence_factory.build(name="good symbols", sequence="ACGT-NRYSWKMDHBV?"),
        sequence_factory.build(name="bad symbols", sequence="L"),
        sequence_factory.build(name="empty seq", sequence=""),
    }


class TestNonnucleotideMatcher:
    def test_affected_sequences(self, all_sequences):
        matcher = Nonnucleotide()
        assert {seq.name: matcher.is_affected(seq) for seq in all_sequences} == {
            "good symbols": True,
            "bad symbols": False,
            "empty seq": True,
        }
