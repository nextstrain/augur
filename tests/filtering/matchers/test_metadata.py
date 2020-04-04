import argparse

from augur.filtering.matchers.metadata import Metadata

import pytest


@pytest.fixture
def all_sequences(sequence_factory):
    return {
        sequence_factory.build(metadata={"color": "red"}),
        sequence_factory.build(metadata={"color": "blue"}),
        sequence_factory.build(metadata={"color": "green"}),
    }


class TestMetadataMatcher:
    @pytest.mark.parametrize(
        "matcher_args, expected_conditions",
        [
            ("color=red", [("color", "=", "red")]),
            ("color!=red", [("color", "!=", "red")]),
            ("color=red,size=med", [("color", "=", "red"), ("size", "=", "med")]),
        ],
    )
    def test_build(self, matcher_args, expected_conditions):
        matcher = Metadata.build(matcher_args)

        assert isinstance(matcher, Metadata)
        assert matcher.conditions == expected_conditions

    @pytest.mark.parametrize(
        "matcher_args, expected_exception",
        [
            ("color=", AttributeError),
            ("color==red", AttributeError),
            ("color", AttributeError),
        ],
    )
    def test_build_bad_arg_string(self, mocker, matcher_args, expected_exception):
        with pytest.raises(expected_exception):
            Metadata.build(matcher_args)

    @pytest.mark.parametrize(
        "matcher_args, expected_matches",
        [
            ("color=red", {"red"}),
            ("color!=red", {"blue", "green"}),
            ("color!=red,color!=blue", {"red", "blue", "green"}),  # yes, an OR of the != conditions covers the entire set
            ("color=red,color=blue", {"red", "blue"}),
            ("power=high", {"red", "blue", "green"}),
        ],
    )
    def test_is_affected(self, all_sequences, matcher_args, expected_matches):
        matcher = Metadata.build(matcher_args)
        assert (
            set(
                [
                    seq.metadata["color"]
                    for seq in all_sequences
                    if matcher.is_affected(seq)
                ]
            )
            == expected_matches
        )
