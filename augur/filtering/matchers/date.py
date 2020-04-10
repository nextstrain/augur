import datetime
import re

from augur.utils import ambiguous_date_to_date_range
from augur.filtering.matchers.base_matcher import BaseMatcher


class Date(BaseMatcher):
    FILTER_OPERATION = "exclude"

    def __init__(self, *, min_date=None, max_date=None):
        if min_date:
            try:
                self.min_date = datetime.date.fromisoformat(min_date)
            except ValueError:
                # TODO is this the legacy behavior? seems weird.
                self.min_date = datetime.date(year=int(min_date), month=1, day=1)
        else:
            self.min_date = None
        if max_date:
            try:
                self.max_date = datetime.date.fromisoformat(max_date)
            except ValueError:
                self.max_date = datetime.date(year=int(max_date), month=1, day=1)
        else:
            self.max_date = None

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

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--min-date", type=str, help="minimal cutoff for numerical date"
        )
        parser.add_argument(
            "--max-date", type=str, help="maximal cutoff for numerical date"
        )
