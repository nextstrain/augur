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
    def test_build(self, mocker):
        mock_read_names = mocker.patch(
            "augur.filtering.matchers.name.read_names",
            mocker.MagicMock(return_value=["a", "b"]),
        )

        matcher_obj = Name.build("file=include.txt")

        mock_read_names.assert_called_once_with("include.txt")
        assert matcher_obj.names == {"a", "b"}

    def test_read_names(self, mocker):
        mocker.patch("builtins.open", mocker.mock_open(read_data="a\nb #comment\n\nc"))

        assert Name.build("file=_").names == {"a", "b", "c"}

    def test_read_names_unsupported_scheme(self, mocker):
        mocker.patch("builtins.open", mocker.mock_open(read_data="a\nb #comment\n\nc"))

        with pytest.raises(Exception):
            assert Name.build("list=a,b.c")

    def test_read_names_empty_file(self, mocker):
        mocker.patch("builtins.open", mocker.mock_open(read_data=""))

        assert Name.build("file=_").names == set()

    def test_read_names_no_such_file(self, mocker):
        with pytest.raises(FileNotFoundError):
            assert Name.build("file=does_not_exist.txt").names == set()

    def test_is_affected_names_in_sequences(self, all_sequences):
        matcher_obj = Name(names=["b", "c"])
        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": True,
            "c": True,
            "d": False,
            "e": False,
        }

    def test_is_affected_names_not_in_sequences(self, all_sequences):
        matcher_obj = Name(names=["g"])
        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": False,
            "c": False,
            "d": False,
            "e": False,
        }

    def test_is_affected_empty_names(self, all_sequences):
        matcher_obj = Name(names=[])
        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": False,
            "c": False,
            "d": False,
            "e": False,
        }

    def test_is_affected_none_names(self, all_sequences):
        matcher_obj = Name(names=[])
        assert {seq.name: matcher_obj.is_affected(seq) for seq in all_sequences} == {
            "a": False,
            "b": False,
            "c": False,
            "d": False,
            "e": False,
        }
