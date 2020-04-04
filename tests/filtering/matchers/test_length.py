import argparse

from augur.filtering.matchers.length import Length

import pytest


@pytest.fixture
def all_sequences(sequence_factory):
    return {
        sequence_factory.build(name="ten", sequence="aaaataaaat"),
        sequence_factory.build(name="five", sequence="aaaat"),
        sequence_factory.build(name="zero", sequence=""),
    }


class TestLengthMatcher:
    def test_build(self, mocker):
        matcher = Length.build("min=5")

        assert isinstance(matcher, Length)
        assert matcher.min_length == 5

    def test_is_affected(self, all_sequences):
        matcher = Length(min_length=5)
        assert {seq.name: matcher.is_affected(seq) for seq in all_sequences} == {
            "ten": False,
            "five": False,
            "zero": True,
        }

    def test_is_affected_none_min_length(self, all_sequences):
        with pytest.raises(TypeError):
            matcher = Length(min_length=None)
            matcher.is_affected(all_sequences[0])
