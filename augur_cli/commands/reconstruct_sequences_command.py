"""TODO reconstruct_sequences command"""

from augur_cli.commands.base_command import BaseCommand
from augur import reconstruct_sequences


class ReconstructSequencesCommand(BaseCommand):
    def run(self):
        reconstruct_sequences.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        reconstruct_sequences.register_arguments(parser)
