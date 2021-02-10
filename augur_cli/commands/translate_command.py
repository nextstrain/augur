"""TODO translate command"""

from augur_cli.commands.base_command import BaseCommand
from augur import translate


class TranslateCommand(BaseCommand):
    def run(self):
        translate.run(self.args)

    @classmethod
    def register_arguments(cls, parser):
        translate.register_arguments(parser)
