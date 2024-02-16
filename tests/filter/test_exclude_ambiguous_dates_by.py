# This file contains functional tests that would normally be written as
# Cram-style tests. However, pytest is nice here since it is easy to use with
# parameterized inputs/outputs (not straightforward to set up for Cram tests¹).
# ¹ https://github.com/nextstrain/augur/pull/1183#discussion_r1142687476

import pytest

from augur.errors import AugurError
from augur.filter._run import run

from . import parse_args, write_metadata


@pytest.mark.parametrize(
    "date, ambiguity",
    [
        # Test complete date strings with ambiguous values.
        ("2019-0X-0X", "any"),
        ("2019-XX-XX", "month"),
        ("2019-XX-XX", "day"),
        ("2019-03-XX", "day"),
        ("201X-XX-XX", "year"),
        ("201X-XX-XX", "month"),
        ("201X-XX-XX", "day"),

        # Test incomplete date strings with ambiguous values.
        ("2019", "month"),
        ("2019", "day"),
        ("2019", "any"),
        ("201X", "year"),
        ("201X", "month"),
        ("201X", "day"),
        ("201X", "any"),
    ],
)
def test_date_is_dropped(tmpdir, date, ambiguity):
    metadata = write_metadata(tmpdir, (("strain","date"),
                                       ("SEQ1"  , date)))
    args = parse_args(f'--metadata {metadata} --exclude-ambiguous-dates-by {ambiguity}')
    with pytest.raises(AugurError, match="All samples have been dropped"):
        run(args)

@pytest.mark.parametrize(
    "date, ambiguity",
    [
        # Test complete date strings without the specified level of ambiguity.
        ("2019-09-03", "any"),
        ("2019-03-XX", "month"),
        ("2019-09-03", "day"),
        ("2019-XX-XX", "year"),

        # Test incomplete date strings without the specified level of ambiguity.
        ("2019", "year"),
    ],
)
def test_date_is_not_dropped(tmpdir, date, ambiguity):
    metadata = write_metadata(tmpdir, (("strain","date"),
                                       ("SEQ1"  , date)))
    args = parse_args(f'--metadata {metadata} --exclude-ambiguous-dates-by {ambiguity}')
    run(args)
