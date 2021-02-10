"""TODO version command"""

from augur_cli.commands.base_command import BaseCommand
from augur import version


class VersionCommand(BaseCommand):
    def run(self):
        version.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        version.register_arguments(parser)
