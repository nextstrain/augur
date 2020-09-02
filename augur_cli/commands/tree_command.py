"""TODO tree command"""

from augur_cli.commands.base_command import BaseCommand
from augur import tree


class TreeCommand(BaseCommand):
    def run(self):
        tree.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        tree.register_arguments(parser)
