from augur.filtering.matchers.base_matcher import BaseMatcher


def read_names(filename):
    with open(filename, "r") as file:
        return [
            line
            for line in (strip_comment(line) for line in file.readlines())
            if line != ""
        ]


def strip_comment(line):
    return line.partition("#")[0].strip()


class Name(BaseMatcher):
    def __init__(self, *, names):
        self.names = set(names)

    def is_affected(self, sequence):
        return sequence.name in self.names

    @classmethod
    def build(cls, matcher_args):
        scheme, _, location = matcher_args.partition("=")

        if scheme == "file":
            return cls(names=read_names(location))
        else:
            raise

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--exclude",
            type=str,
            help="file with list of strains that are to be excluded",
        )
        parser.add_argument(
            "--include",
            type=str,
            help="file with list of strains that are to be included regardless of priorities or subsampling",
        )
