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
