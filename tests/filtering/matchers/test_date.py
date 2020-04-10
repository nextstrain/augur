import argparse
import datetime

from augur.filtering.matchers.date import Date

import pytest


@pytest.fixture
def all_sequences(sequence_factory):
    return {
        sequence_factory.build(metadata={"date": "1991-3-1"}),
        sequence_factory.build(metadata={"date": "2000-3-10"}),
        sequence_factory.build(metadata={"date": "2000-3-XX"}),
        sequence_factory.build(metadata={"date": "2009-3-1"}),
        sequence_factory.build(metadata={}),
    }


class TestDateMatcher:
    @pytest.mark.parametrize(
        "min_arg, max_arg, expected_min, expected_max",
        [
            ("2000-01-01", None, datetime.date(year=2000, month=1, day=1), None),
            (None, "2000-02-02", None, datetime.date(year=2000, month=2, day=2)),
            (
                "2000-03-03",
                "2000-04-04",
                datetime.date(year=2000, month=3, day=3),
                datetime.date(year=2000, month=4, day=4),
            ),
            (
                "2003",
                "2004",
                datetime.date(year=2003, month=1, day=1),
                datetime.date(year=2004, month=1, day=1),
            ),
        ],
    )
    def test_init(self, mocker, min_arg, max_arg, expected_min, expected_max):
        matcher = Date(min_date=min_arg, max_date=max_arg)

        assert isinstance(matcher, Date)
        assert matcher.min_date == expected_min
        assert matcher.max_date == expected_max

    def test_is_affected_min_only(self, all_sequences):
        matcher = Date(min_date="2000-03-22")

        assert {seq.date: matcher.is_affected(seq) for seq in all_sequences} == {
            "1991-3-1": True,
            "2000-3-10": True,
            "2000-3-XX": False,
            "2009-3-1": False,
            None: False,
        }

    def test_is_affected_max_only(self, all_sequences):
        matcher = Date(max_date="2000-03-05")

        assert {seq.date: matcher.is_affected(seq) for seq in all_sequences} == {
            "1991-3-1": False,
            "2000-3-10": True,
            "2000-3-XX": False,
            "2009-3-1": True,
            None: False,
        }

    def test_is_affected_min_and_max(self, all_sequences):
        matcher = Date(min_date="2000-03-20", max_date="2000-03-22")

        assert {seq.date: matcher.is_affected(seq) for seq in all_sequences} == {
            "1991-3-1": True,
            "2000-3-10": True,
            "2000-3-XX": False,
            "2009-3-1": True,
            None: False,
        }

    def test_is_affected_both_none(self, all_sequences):
        matcher = Date(min_date=None, max_date=None)

        assert {seq.date: matcher.is_affected(seq) for seq in all_sequences} == {
            "1991-3-1": False,
            "2000-3-10": False,
            "2000-3-XX": False,
            "2009-3-1": False,
            None: False,
        }
