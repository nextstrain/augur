"""TODO mask command"""

from augur_cli.commands.base_command import BaseCommand
from augur import align


class AlignCommand(BaseCommand):
    def run(self):
        align.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        align.register_arguments(parser)
