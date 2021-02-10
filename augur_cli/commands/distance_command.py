"""TODO distance command"""

from augur_cli.commands.base_command import BaseCommand
from augur import distance


class DistanceCommand(BaseCommand):
    def run(self):
        distance.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        distance.register_arguments(parser)
