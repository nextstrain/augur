"""TODO mask command"""

from augur_cli.commands.base_command import BaseCommand
from augur import mask


class MaskCommand(BaseCommand):
    def run(self):
        mask.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        mask.register_arguments(parser)
