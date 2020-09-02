"""TODO traits command"""

from augur_cli.commands.base_command import BaseCommand
from augur import traits


class TraitsCommand(BaseCommand):
    def run(self):
        traits.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        traits.register_arguments(parser)
