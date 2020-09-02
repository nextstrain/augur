import argparse
import os
import sys
from types import SimpleNamespace

from augur_cli.commands import COMMAND_CLASSES

recursion_limit = os.environ.get("AUGUR_RECURSION_LIMIT")
if recursion_limit:
    sys.setrecursionlimit(int(recursion_limit))


def make_parser():
    parser = argparse.ArgumentParser(
        prog="augur",
        description="Augur: A bioinformatics toolkit for phylogenetic analysis.",
    )

    add_default_command(parser)
    add_version_alias(parser)

    subparsers = parser.add_subparsers()
    for command_class in COMMAND_CLASSES:
        command_name = command_class.command_name()
        command_desc = sys.modules[command_class.__module__].__doc__.strip()

        subparser = subparsers.add_parser(
            command_name, help=command_desc.splitlines()[0], description=command_desc
        )

        subparser.set_defaults(__command_class__=command_class)
        subparser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

        command_class.register_arguments(subparser)

    return parser


def main():
    args = make_parser().parse_args(sys.argv[1:])
    try:
        return args.__command_class__(args).run()
    except RecursionError:
        print(
            f"FATAL: Maximum recursion depth reached. You can set the env variable AUGUR_RECURSION_LIMIT to adjust this (current limit: {sys.getrecursionlimit()})"
        )
        sys.exit(2)


def add_default_command(parser):
    """
    Sets the default command to run when none is provided.
    """

    class default_command:
        def run(args):
            parser.print_help()
            return 2

    parser.set_defaults(__command_class__=default_command)


def add_version_alias(parser):
    """
    Add --version as a (hidden) alias for the version command.

    It's not uncommon to blindly run a command with --version as the sole
    argument, so its useful to make that Just Work.
    """

    class run_version_command(argparse.Action):
        def __call__(self, *args, **kwargs):
            opts = SimpleNamespace()
            sys.exit(version.run(opts))

    return parser.add_argument(
        "--version", nargs=0, help=argparse.SUPPRESS, action=run_version_command
    )


def command_name(command_module):
    """
    Returns a short name for a command module.
    TODO move to BaseCommand class
    """

    return command_module.__name__.split(".")[-1].replace("_", "-")


if __name__ == "__main__":
    exit(main())
