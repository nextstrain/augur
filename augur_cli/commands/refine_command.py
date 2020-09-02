"""TODO refine command"""

from augur_cli.commands.base_command import BaseCommand
from augur import refine


class RefineCommand(BaseCommand):
    def run(self):
        refine.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        refine.register_arguments(parser)
