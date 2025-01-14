# Originally based on tests/help.py from the Nextstrain CLI project.¹
#
# ¹ <https://github.com/nextstrain/cli/blob/4a00d7100eff811eab6df34db73c7f6d4196e22b/tests/help.py>
import pytest

from augur import make_parser
from augur.argparse_ import walk_commands, ExtendOverwriteDefault


# Walking commands is slow, so do it only once and use it for all tests in this
# file (though currently n=1).
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
