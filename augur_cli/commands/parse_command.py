"""TODO parse command"""

from augur_cli.commands.base_command import BaseCommand
from augur import parse


class ParseCommand(BaseCommand):
    def run(self):
        parse.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        parse.register_arguments(parser)
