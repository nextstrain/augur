import abc


class BaseMatcher(abc.ABC):
    """
    def filter(self, sequences, all_sequences):
        if self.filter_operation() == "include":
            self.affected_sequences_set = set(self.affected_sequences(all_sequences))

            return set(sequences).union(self.affected_sequences_set)
        elif self.filter_operation() == "exclude":
            self.affected_sequences_set = set(self.affected_sequences(sequences))

            return set(sequences).difference(self.affected_sequences_set)
        else:
            raise Exception(f"Unknown filter operation {self.filter_operation()}")
    """

    @abc.abstractmethod
    def is_affected(self, sequence):
        pass

    """
    @staticmethod
    @abc.abstractmethod
    def is_enabled(args):
        pass
    """
