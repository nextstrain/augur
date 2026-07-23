"""
Custom helpers for the argparse standard library.
"""
import os
import sys
import re
import textwrap
import configargparse
from argparse import Action, ArgumentDefaultsHelpFormatter, _ArgumentGroup, _SubParsersAction
from itertools import chain
from collections import OrderedDict
from textwrap import dedent, indent as indent_text
from typing import Iterable, Optional, Tuple, Union
from configargparse import YAMLConfigFileParser
from .rst import rst_to_text
from .types import ValidationMode


class AugurYAMLConfigFileParser(YAMLConfigFileParser):
    """
    Custom YAML parser that maps config keys with underscores to dashes.
    This allows users to specify options like `use_fft: true` in the YAML
    config, while mapping internally to the `--use-fft` argparse flag.
    """
    def parse(self, stream):
        parsed = super().parse(stream)
        new_parsed = OrderedDict()
        for k, v in parsed.items():
            new_parsed[k.replace('_', '-')] = v
        return new_parsed


def _is_option_in_args(option_string, args_list):
    for arg in args_list:
        if arg == option_string or arg.startswith(option_string + "="):
            return True
        if option_string.startswith("-") and not option_string.startswith("--") and len(option_string) == 2 and arg.startswith(option_string):
            return True
    return False


class CustomArgumentParser(configargparse.ArgumentParser):
    """
    Subclass of configargparse.ArgumentParser that prevents command-line arguments
    from overriding config file values. If an argument is provided both in a config file
    and on the command line, an error is raised.
    """
    def parse_known_args(self, args=None, namespace=None, **kwargs):
        parsed_args, argv = super().parse_known_args(args=args, namespace=namespace, **kwargs)

        config_file_path = getattr(parsed_args, "config", None)
        if config_file_path:
            try:
                with open(config_file_path) as fp:
                    config_keys = set(self._config_file_parser.parse(fp).keys())
            except Exception:
                config_keys = set()

            args_list = sys.argv[1:] if args is None else args
            cli_dests = set()
            for action in self._actions:
                if getattr(action, "is_config_file_arg", False):
                    continue
                for opt in action.option_strings:
                    if _is_option_in_args(opt, args_list):
                        cli_dests.add(opt.lstrip("-"))
                        break

            conflicts = sorted(config_keys & cli_dests)
            if conflicts:
                self.error(
                    f"The following option(s) were specified both in the config file and on the CLI: "
                    f"{', '.join(conflicts)}"
                )

        return parsed_args, argv

    def format_help(self):
        help_text = super().format_help()
        pattern = r"In general, command-line values\s+override config file values which override defaults\."
        replacement = "Specifying an option on both the command line and in a config file will result in an error."
        help_text = re.sub(pattern, replacement, help_text)

        if "Args that start with" in help_text:
            main_help, footer = help_text.rsplit("Args that start with", 1)
            footer = "Args that start with" + footer
            text_width = max(self._get_formatter()._width, 11)
            footer_clean = " ".join(footer.split())
            refilled = textwrap.fill(footer_clean, text_width)
            help_text = main_help + refilled + "\n"

        return help_text


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


def add_subparser(parent_subparsers: _SubParsersAction, *args, **kwargs) -> CustomArgumentParser:
    """
    Add a subparser to a parent subparser.
    """
    # Use the same formatting class for every command for consistency.
    kwargs.setdefault("formatter_class", HelpFormatter)

    return parent_subparsers.add_parser(*args, **kwargs)


SUBPARSER_ATTRIBUTE = '__subparser__'

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
        
        # Ensure subparsers use the custom YAML parser
        subparser._config_file_parser = AugurYAMLConfigFileParser()

        subparser.set_defaults(**{
            # Add default attribute for command module
            command_attribute: command,

            # Add a reference to the subparser
            SUBPARSER_ATTRIBUTE: subparser,
        })

        # Use the same formatting class for every command for consistency.
        # Set here to avoid repeating it in every command's register_parser().
        subparser.formatter_class = HelpFormatter

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


def add_validation_arguments(parser: Union[CustomArgumentParser, _ArgumentGroup]):
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
def walk_commands(parser: CustomArgumentParser, command: Optional[Tuple[str, ...]] = None) -> Iterable[Tuple[Tuple[str, ...], CustomArgumentParser]]:
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
