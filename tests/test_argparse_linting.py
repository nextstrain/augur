# Originally based on tests/help.py from the Nextstrain CLI project.¹
#
# ¹ <https://github.com/nextstrain/cli/blob/4a00d7100eff811eab6df34db73c7f6d4196e22b/tests/help.py>
import pytest
import sys

from augur import make_parser
from augur.argparse_ import walk_commands, ExtendOverwriteDefault
from subprocess import run


# Walking commands is slow, so do it only once and use it for all tests in this
# file.
commands = list(walk_commands(make_parser()))


def help_texts():
    for command, parser in commands:
        for attribute in ("description", "epilog"):
            yield command, attribute, getattr(parser, attribute)

        for group in parser._action_groups:
            yield command, f"{group.title} description", group.description

        for action in parser._actions:
            yield command, "/".join(action.option_strings) or action.dest, action.help


# Ensure we always use ExtendOverwriteDefault for options that take a variable
# number of arguments.  See <https://github.com/nextstrain/augur/pull/1709>.
@pytest.mark.parametrize("action", [
    pytest.param(action, id = " ".join(command) + " " + "/".join(action.option_strings))
        for command, parser in commands
        for action in parser._actions
         if action.nargs in {"+", "*"}
])
def test_ExtendOverwriteDefault(action):
    assert isinstance(action, ExtendOverwriteDefault)


@pytest.mark.parametrize("command", [command for command, parser in commands], ids = lambda command: " ".join(command))
def test_help(command):
    args = [sys.executable, "-m", "augur", *command[1:], "--help"]

    # Check the exit status ourselves for nicer test output on failure
    result = run(args)
    assert result.returncode == 0, f"{args} exited with error"

    # TODO: Test with AUGUR_RST_STRICT once all Augur help strings are valid strict rST.


@pytest.mark.parametrize("command, location, text", [
    pytest.param(command, location, text, id = " ".join(command) + " " + location)
        for command, location, text in help_texts()
])
def test_multiline_help_dedented(command, location, text):
    assert not _has_extra_indent(text), \
        f"{' '.join(command)} {location} has extra indentation"


def _has_extra_indent(text):
    # Ignore non-strings.
    if not isinstance(text, str):
        return False

    # Ignore blank lines, since leading/trailing blank lines are common in
    # triple-quoted strings and do not affect indentation.
    lines = [line for line in text.splitlines() if line.strip()]

    if len(lines) < 2:
        return False

    indents = [
        len(line) - len(line.lstrip())
        for line in lines
    ]

    return (
        # Fail if all lines after the first are indented.
        all(indent >= 1 for indent in indents[1:])
    )
