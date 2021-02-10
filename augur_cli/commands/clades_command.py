"""TODO clades command"""

from augur_cli.commands.base_command import BaseCommand
from augur import clades


class CladesCommand(BaseCommand):
    def run(self):
        clades.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        clades.register_arguments(parser)
