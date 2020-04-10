import re

from augur.filtering.matchers.base_matcher import BaseMatcher


class Length(BaseMatcher):
    def __init__(self, *, min_length):
        self.min_length = int(min_length)

    def is_affected(self, sequence):
        return sequence.length < self.min_length

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--min-length", type=int, help="minimal length of the sequences"
        )
