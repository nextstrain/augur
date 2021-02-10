"""TODO frequencies command"""

from augur_cli.commands.base_command import BaseCommand
from augur import frequencies


class FrequenciesCommand(BaseCommand):
    def run(self):
        frequencies.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        frequencies.register_arguments(parser)
