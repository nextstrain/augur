"""TODO ancestral command"""

from augur_cli.commands.base_command import BaseCommand
from augur import ancestral


class AncestralCommand(BaseCommand):
    def run(self):
        ancestral.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        ancestral.register_arguments(parser)
