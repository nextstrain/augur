"""TODO export command"""

from augur_cli.commands.base_command import BaseCommand
from augur import export


class ExportCommand(BaseCommand):
    def run(self):
        export.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        export.register_arguments(parser)
