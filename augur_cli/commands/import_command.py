"""
Import analyses into augur pipeline from other systems
"""

from augur_cli.commands.base_command import BaseCommand
from augur import import_beast


class ImportCommand(BaseCommand):
    def run(self):
        if "beast" in self.args:
            return import_beast.run_beast(self.args)

    @classmethod
    def register_arguments(cls, parser):
        metavar_msg = "Import analyses into augur pipeline from other systems"
        subparsers = parser.add_subparsers(title="TYPE", metavar=metavar_msg)
        subparsers.required = True
        import_beast.register_arguments_beast(subparsers)
