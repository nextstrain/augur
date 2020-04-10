import re

from augur.filtering.matchers.base_matcher import BaseMatcher


class Metadata(BaseMatcher):
    def __init__(self, *, clause):
        self.conditions = []
        for condition in clause.split(","):
            matches = re.search(r"(\w+)(=|!=)(\w+)", condition)
            if matches is None:
                raise AttributeError(
                    f"Malformed condition`{condition}`. "
                    "Must be of form `property=value` or `property!=value`"
                )
            self.conditions.append(matches.groups())

    def is_affected(self, sequence):
        return any(
            [
                self.sequence_matches_condition(sequence, key, value, comparator)
                for key, comparator, value in self.conditions
            ]
        )

    @staticmethod
    def sequence_matches_condition(sequence, key, value, comparator):
        if key not in sequence.metadata:
            # sequences that lack the specified metadata column are always excluded
            # TODO is that right? regardless, it should be documented in the argument help
            return True

        sequence_value = sequence.metadata[key]
        if comparator == "=":
            return sequence_value.lower() == value.lower()

        if comparator == "!=":
            return sequence_value.lower() != value.lower()

        # Should be impossible
        raise

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            "--exclude-where",
            nargs="+",
            help="Exclude samples matching these conditions. Ex: `host=rat` or "
            "`host!=rat`. Multiple values are processed as OR (matching any of "
            "those specified will be excluded), not AND",
        )
        parser.add_argument(
            "--include-where",
            nargs="+",
            help="Include samples with these values. ex: host=rat. Multiple values "
            "are processed as OR (having any of those specified will be included), "
            "not AND. This rule is applied last and ensures any sequences matching "
            "these rules will be included.",
        )
