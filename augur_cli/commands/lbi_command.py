"""TODO lbi command"""

from augur_cli.commands.base_command import BaseCommand
from augur import lbi


class LbiCommand(BaseCommand):
    def run(self):
        lbi.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        lbi.register_arguments(parser)
