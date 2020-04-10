from augur.filtering.matchers.base_matcher import BaseMatcher


class Nonnucleotide(BaseMatcher):
    def is_affected(self, sequence):
        return not sequence.has_only_nucleotide_symbols

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--non-nucleotide",
            action="store_true",
            help="exclude sequences that contain illegal characters",
        )
