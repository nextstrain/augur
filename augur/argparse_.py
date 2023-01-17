"""
Custom helpers for the argparse standard library.
"""
from argparse import Action, ArgumentDefaultsHelpFormatter

from augur.utils import first_line


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

    commands: list[ModuleType]
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
            subparser.description = first_line(command.__doc__)

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


# XXX TODO: Drop this when we drop support for 3.7 and use builtin "extend"
# action.
class ExtendAction(Action):
    """
    Loosely backported version of builtin "extend" action
    (argparse._ExtendAction) in CPython v3.10.5-172-gf118661a18.  Some of the
    ceremony we don't need has been ditched.

    >>> import argparse
    >>> p = argparse.ArgumentParser()
    >>> a = p.add_argument("-x", action=ExtendAction, nargs="+")
    >>> p.parse_args(["-x", "a", "b", "-x", "c", "-x", "d"]).x
    ['a', 'b', 'c', 'd']
    """
    def __call__(self, parser, namespace, values, option_string=None):
        current = getattr(namespace, self.dest, None) or []
        setattr(namespace, self.dest, [*current, *values])
