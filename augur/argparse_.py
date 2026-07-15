"""
Custom helpers for the argparse standard library.
"""
import os
from argparse import Action, ArgumentDefaultsHelpFormatter, ArgumentParser, _ArgumentGroup, _SubParsersAction
from itertools import chain
from textwrap import dedent, indent as indent_text
from typing import Iterable, Optional, Tuple, Union
from .rst import rst_to_text
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


class HelpFormatter(ArgumentDefaultsHelpFormatter):
    def __init__(self, prog, indent_increment = 2, max_help_position = 24, width = None):
        # Ignore terminal size, unlike standard argparse, as the readability of
        # paragraphs of text suffers at wide widths.  Instead, default to 78
        # columns (80 wide - 2 column gutter), but let that be overridden by an
        # explicit COLUMNS variable or reduced by a smaller actual terminal.
        if width is None:
            try:
                width = int(os.environ["COLUMNS"])
            except (KeyError, ValueError):
                try:
                    width = min(os.get_terminal_size().columns, 80) - 2
                except (AttributeError, OSError):
                    width = 80 - 2

        super().__init__(prog, indent_increment, max_help_position, width)

    # Based on argparse.RawDescriptionHelpFormatter's implementation
    def _fill_text(self, text, width, indent):
        return indent_text(rst_to_text(text, width), prefix=indent)

    # Based on argparse.RawTextHelpFormatter's implementation
    def _split_lines(self, text, width):
        # Render to rST here so rST gets control over wrapping/line breaks.
        return rst_to_text(text, width).splitlines()

    # Emit blank lines between arguments for better readability.  It might seem
    # simpler to override _join_parts() instead, but that's called from two
    # places and we only want to join on newlines for one of those places.
    def add_arguments(self, actions):
        for i, action in enumerate(actions):
            if i != 0:
                # Use " \n" to avoid a completely empty line (e.g. "\n") for
                # the sake ShowBriefHelp.truncate_help()'s heuristic.
                self._add_item(str, [" \n"])

            self.add_argument(action)


def add_default_command(parser):
    """
    Sets the default command to run when none is provided.
    """
    class default_command():
        def run(args):
            parser.print_help()
            return 2

    parser.set_defaults(__command__ = default_command)


def add_subparser(parent_subparsers: _SubParsersAction, *args, **kwargs) -> ArgumentParser:
    """
    Add a subparser to a parent subparser.
    """
    # Use the same formatting class for every command for consistency.
    kwargs.setdefault("formatter_class", HelpFormatter)

    return parent_subparsers.add_parser(*args, **kwargs)


def register_commands(parser: ArgumentParser, commands, command_attribute='__command__'):
    """
    Add subparsers for each command module.

    Parameters
    ----------
    parser
        ArgumentParser object.
    commands: list[types.ModuleType]
        A list of modules that are commands that require their own subparser.
        Each module is required to have a `register_parser` function to add its own
        subparser and arguments.

    command_attribute: str, optional
        Optional attribute name for the commands. The default is `__command__`,
        which allows top level augur to run commands directly via `args.__command__.run()`.
    """
    subparsers = parser.add_subparsers(dest="subcommand", required=False)

    for command in commands:
        # Allow each command to register its own subparser
        subparser = command.register_parser(subparsers)

        # Add default attribute for command module
        if command_attribute:
            subparser.set_defaults(**{command_attribute: command})

        # Ensure all subparsers format like the top-level parser
        subparser.formatter_class = parser.formatter_class

        if not subparser.description and command.__doc__:
            subparser.description = command.__doc__

        # If a command doesn't have its own run() function, then print its help when called.
        if not getattr(command, "run", None):
            add_default_command(subparser)

        # Recursively register any subcommands
        if getattr(subparser, "subcommands", None):
            register_commands(subparser, subparser.subcommands)


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
        help=dedent("""\
            Control if optional validation checks are performed and what
            happens if they fail.

            'error' and 'warn' modes perform validation and emit messages about
            failed validation checks.  'error' mode causes a non-zero exit
            status if any validation checks failed, while 'warn' does not.

            'skip' mode performs no validation.

            Note that some validation checks are non-optional and as such are
            not affected by this setting.
        """))
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
