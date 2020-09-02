"""
Validate files related to augur consumption or export.
"""

import sys

from augur_cli.commands.base_command import BaseCommand
import augur


class ValidateCommand(BaseCommand):
    def run(self):
        pass

    def auspice_config_v2(self, args):
        self.validate_file(args.fname, validator=augur.validate.auspice_config_v2)

    def export_v2(self, args):
        self.validate_file(args.fname, validator=augur.validate.export_v2)

    def export_v1(self, args):
        self.validate_file(
            args.meta_fname, args.tree_fname, validator=augur.validate.export_v1
        )

    def validate_file(self, fname, *, validator):
        print(f"Validating schema of {fname}...", end="")

        try:
            validator(fname)
        except augur.validate.JsonValidationError as e:
            for error in e.errors:
                print(error, file=sys.stderr)

        print("Passed")

    @classmethod
    def register_arguments(cls, parser):
        subparsers = parser.add_subparsers(
            dest="subcommand", help="File(s) to validate"
        )

        subparser = subparsers.add_parser(
            "export-v2", help="validate JSON intended for auspice v2"
        )
        subparser.add_argument(
            "fname", metavar="FILE", help="exported (main) v2 auspice JSON"
        )

        subparser = subparsers.add_parser(
            "export-v1", help="validate tree+meta JSONs intended for auspice v1"
        )
        subparser.add_argument(
            "meta_fname", metavar="FILE", help="exported (v1) meta JSON"
        )
        subparser.add_argument(
            "tree_fname", metavar="FILE", help="exported (v1) tree JSON"
        )

        subparser = subparsers.add_parser(
            "auspice-config-v2",
            help="validate auspice config intended for `augur export v2`",
        )
        subparser.add_argument("fname", metavar="JSON", help="auspice config JSON")
