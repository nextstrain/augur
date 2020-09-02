"""TODO sequence_traits command"""

from augur_cli.commands.base_command import BaseCommand
from augur import sequence_traits


class SequenceTraitsCommand(BaseCommand):
    def run(self):
        sequence_traits.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        sequence_traits.register_arguments(parser)
