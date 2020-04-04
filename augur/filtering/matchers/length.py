import re

from augur.filtering.matchers.base_matcher import BaseMatcher


class Length(BaseMatcher):
    def __init__(self, *, min_length):
        self.min_length = int(min_length)

    def is_affected(self, sequence):
        return sequence.length < self.min_length

    @classmethod
    def build(cls, matcher_args):
        arg_matches = re.search(r"min=([^,]*)", matcher_args)
        min_length = arg_matches.groups()[0] if arg_matches else None

        return cls(min_length=min_length)

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--min-length", type=int, help="minimal length of the sequences"
        )
