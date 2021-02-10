"""TODO titers command"""

from augur_cli.commands.base_command import BaseCommand
from augur import titers


class TitersCommand(BaseCommand):
    def run(self):
        titers.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        titers.register_arguments(parser)
