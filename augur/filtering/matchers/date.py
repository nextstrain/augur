import datetime
import re

from augur.utils import ambiguous_date_to_date_range
from augur.filtering.matchers.base_matcher import BaseMatcher


class Date(BaseMatcher):
    FILTER_OPERATION = "exclude"

    def __init__(self, *, min_date=None, max_date=None):
        self.min_date = min_date
        self.max_date = max_date

    def is_affected(self, sequence):
        # TODO need the ambiguous sample date be ENTIRELY within the constraint range? or is partial overlap enough to inclued the sample?
        if sequence.date is None:
            return False

        sample_date_min, sample_date_max = ambiguous_date_to_date_range(
            sequence.date, "%Y-%m-%d"
        )

        if self.min_date and sample_date_max < self.min_date:
            return True

        if self.max_date and sample_date_min > self.max_date:
            return True

        return False

    @classmethod
    def build(cls, matcher_args):
        # TODO naive and clumsy
        arg_matches = re.search(r"min=([^,]*)", matcher_args)
        min_date_arg = arg_matches.groups()[0] if arg_matches else None

        arg_matches = re.search(r"max=([^,]*)", matcher_args)
        max_date_arg = arg_matches.groups()[0] if arg_matches else None

        # TODO fromisoformat is a behavior change! might as well introduce some predictability to the input. backwards compatibility?
        if min_date_arg:
            try:
                min_date = datetime.date.fromisoformat(min_date_arg)
            except ValueError:
                # TODO is this the legacy behavior? seems weird.
                min_date = datetime.date(year=int(min_date_arg), month=1, day=1)
        else:
            min_date = None
        if max_date_arg:
            try:
                max_date = datetime.date.fromisoformat(max_date_arg)
            except ValueError:
                max_date = datetime.date(year=int(max_date_arg), month=1, day=1)
        else:
            max_date = None
        return Date(min_date=min_date, max_date=max_date)

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--min-date", type=str, help="minimal cutoff for numerical date"
        )
        parser.add_argument(
            "--max-date", type=str, help="maximal cutoff for numerical date"
        )
