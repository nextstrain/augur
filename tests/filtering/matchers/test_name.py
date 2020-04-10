import argparse

from augur.filtering.matchers.name import Name

import pytest


@pytest.fixture
def all_sequences(sequence_factory):
    return {
        sequence_factory.build(name="a"),
        sequence_factory.build(name="b"),
        sequence_factory.build(name="c"),
        sequence_factory.build(name="d"),
        sequence_factory.build(name="e"),
    }


class TestNameMatcher:
    def test_init(self, mocker):
        mocker.patch("builtins.open", mocker.mock_open(read_data="a\nb"))

        matcher_obj = Name(filename="_")

        assert matcher_obj.names == {"a", "b"}

    def test_read_names(self, mocker):
        mocker.patch("builtins.open", mocker.mock_open(read_data="a\nb #comment\n\nc"))

        assert Name(filename="_").names == {"a", "b", "c"}

    def test_read_names_none_file(self, mocker):
        mocker.patch("builtins.open", mocker.mock_open(read_data=""))

        assert Name(filename=None).names == set()

    def test_read_names_no_such_file(self, mocker):
        with pytest.raises(FileNotFoundError):
            assert Name(filename="does_not_exist.txt").names == set()

    def test_is_affected_names_in_sequences(self, mocker, all_sequences):
        mocker.patch("augur.filtering.matchers.name.read_names", lambda _: ["b", "c"])

        matcher_obj = Name(filename="_")

        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": True,
            "c": True,
            "d": False,
            "e": False,
        }

    def test_is_affected_names_not_in_sequences(self, mocker, all_sequences):
        mocker.patch("augur.filtering.matchers.name.read_names", lambda _: ["g"])

        matcher_obj = Name(filename="_")

        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": False,
            "c": False,
            "d": False,
            "e": False,
        }

    def test_is_affected_empty_names(self, mocker, all_sequences):
        mocker.patch("augur.filtering.matchers.name.read_names", lambda _: [])

        matcher_obj = Name(filename="_")

        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": False,
            "c": False,
            "d": False,
            "e": False,
        }
