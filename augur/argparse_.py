"""
Custom helpers for the argparse standard library.
"""
from argparse import Action, ArgumentDefaultsHelpFormatter, ArgumentParser, _ArgumentGroup, _SubParsersAction
from itertools import chain
from typing import Iterable, Optional, Tuple, Union
from .types import ValidationMode


# Include this in an argument help string to suppress the automatic appending
# of the default value by argparse.ArgumentDefaultsHelpFormatter.  This works
# because the automatic appending is conditional on the presence of %(default),
# so we include it but then format it as a zero-length string .0s.  🙃
#
# Another solution would be to add an extra attribute to the argument (the
# argparse.Action instance) and then subclass ArgumentDefaultsHelpFormatter to
# condition on that new attribute, but that seems more brittle.
#
# Copied from the Nextstrain CLI repo
# https://github.com/nextstrain/cli/blob/017c53805e8317951327d24c04184615cc400b09/nextstrain/cli/argparse.py#L13-L21
SKIP_AUTO_DEFAULT_IN_HELP = "%(default).0s"


def add_default_command(parser):
    """
    Sets the default command to run when none is provided.
    """
    class default_command():
        def run(args):
            parser.print_help()
            return 2

    parser.set_defaults(__command__ = default_command)


def add_command_subparsers(subparsers, commands, command_attribute='__command__'):
    """
    Add subparsers for each command module.

    Parameters
    ----------
    subparsers: argparse._SubParsersAction
        The special subparsers action object created by the parent parser
        via `parser.add_subparsers()`.

    commands: list[types.ModuleType]
        A list of modules that are commands that require their own subparser.
        Each module is required to have a `register_parser` function to add its own
        subparser and arguments.

    command_attribute: str, optional
        Optional attribute name for the commands. The default is `__command__`,
        which allows top level augur to run commands directly via `args.__command__.run()`.
    """
    for command in commands:
        # Allow each command to register its own subparser
        subparser = command.register_parser(subparsers)

        # Add default attribute for command module
        if command_attribute:
            subparser.set_defaults(**{command_attribute: command})

        # Use the same formatting class for every command for consistency.
        # Set here to avoid repeating it in every command's register_parser().
        subparser.formatter_class = ArgumentDefaultsHelpFormatter

        if not subparser.description and command.__doc__:
            subparser.description = command.__doc__

        # If a command doesn't have its own run() function, then print its help when called.
        if not getattr(command, "run", None):
            add_default_command(subparser)


class HideAsFalseAction(Action):
    """
    Custom argparse Action that stores False for arguments passed as `--hide*`
    and stores True for all other argument patterns.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, option_string[2:6] != 'hide')


class ExtendOverwriteDefault(Action):
    """
    Similar to the core argparse ``extend`` action, but overwrites the argument
    ``default``, if any, instead of appending to it.

    Thus, the ``default`` value is not included when the option is given and
    may be a non-list value if desired.
    """
    def __call__(self, parser, namespace, value, option_string = None):
        current = getattr(namespace, self.dest, None)

        if current is parser.get_default(self.dest) or current is None:
            current = []

        setattr(namespace, self.dest, [*current, *value])


def add_validation_arguments(parser: Union[ArgumentParser, _ArgumentGroup]):
    """
    Add arguments to configure validation mode of node data JSON files.
    """
    parser.add_argument(
        '--validation-mode',
        dest="validation_mode",
        type=ValidationMode,
        choices=[mode for mode in ValidationMode],
        default=ValidationMode.ERROR,
        help="""
            Control if optional validation checks are performed and what
            happens if they fail.

            'error' and 'warn' modes perform validation and emit messages about
            failed validation checks.  'error' mode causes a non-zero exit
            status if any validation checks failed, while 'warn' does not.

            'skip' mode performs no validation.

            Note that some validation checks are non-optional and as such are
            not affected by this setting.
        """)
    parser.add_argument(
        '--skip-validation',
        dest="validation_mode",
        action="store_const",
        const=ValidationMode.SKIP,
        help="Skip validation of input/output files, equivalent to --validation-mode=skip. Use at your own risk!")


# Originally copied from nextstrain/cli/argparse.py in the Nextstrain CLI project¹.
#
# ¹ <https://github.com/nextstrain/cli/blob/4a00d7100eff811eab6df34db73c7f6d4196e22b/nextstrain/cli/argparse.py#L252-L271>
def walk_commands(parser: ArgumentParser, command: Optional[Tuple[str, ...]] = None) -> Iterable[Tuple[Tuple[str, ...], ArgumentParser]]:
    if command is None:
        command = (parser.prog,)

    yield command, parser

    subparsers = chain.from_iterable(
        action.choices.items()
            for action in parser._actions
             if isinstance(action, _SubParsersAction))

    visited = set()

    for subname, subparser in subparsers:
        if subparser in visited:
            continue

        visited.add(subparser)

        yield from walk_commands(subparser, (*command, subname))
