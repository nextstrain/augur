"""TODO filter command"""

from augur_cli.commands.base_command import BaseCommand
from augur import filter  #TODO builtin conflict


class FilterCommand(BaseCommand):
    def run(self):
        filter.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        filter.register_arguments(parser)
