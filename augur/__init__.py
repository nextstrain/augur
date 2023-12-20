"""
The top-level augur command which dispatches to subcommands.
"""

import argparse
import os
import sys
import importlib
import traceback
from textwrap import dedent
from types import SimpleNamespace
from treetime import TreeTimeError, TreeTimeUnknownError

from .errors import AugurError
from .io.print import print_err
from .argparse_ import add_command_subparsers, add_default_command

DEFAULT_AUGUR_RECURSION_LIMIT = 10000
sys.setrecursionlimit(int(os.environ.get("AUGUR_RECURSION_LIMIT") or DEFAULT_AUGUR_RECURSION_LIMIT))

command_strings = [
    "parse",
    "curate",
    "index",
    "filter",
    "mask",
    "align",
    "tree",
    "refine",
    "ancestral",
    "translate",
    "reconstruct_sequences",
    "clades",
    "traits",
    "sequence_traits",
    "lbi",
    "distance",
    "titers",
    "frequencies",
    "export",
    "validate",
    "version",
    "import_",
    "measurements",
]

COMMANDS = [importlib.import_module('augur.' + c) for c in command_strings]

def make_parser():
    parser = argparse.ArgumentParser(
        prog        = "augur",
        description = "Augur: A bioinformatics toolkit for phylogenetic analysis.")

    add_default_command(parser)
    add_version_alias(parser)

    subparsers = parser.add_subparsers()
    add_command_subparsers(subparsers, COMMANDS)

    return parser


def run(argv):
    args = make_parser().parse_args(argv)
    try:
        return args.__command__.run(args)
    except AugurError as e:
        print_err(f"ERROR: {e}")
        sys.exit(2)
    except RecursionError:
        print_err("FATAL: Maximum recursion depth reached. You can set the env variable AUGUR_RECURSION_LIMIT to adjust this (current limit: {})".format(sys.getrecursionlimit()))
        sys.exit(2)
    except FileNotFoundError as e:
        print_err(f"ERROR: {e.strerror}: '{e.filename}'")
        sys.exit(2)
    except TreeTimeUnknownError as e:
        print_err(dedent("""\
            ERROR from TreeTime: An error occurred in TreeTime (see above). This may be due to an issue with TreeTime or Augur.
            Please report you are calling TreeTime via Augur.
            """))
        sys.exit(2)
    except TreeTimeError as e:
        print_err(f"ERROR: {e}")
        print_err("\n")
        print_err(dedent("""\
            ERROR from TreeTime: This error is most likely due to a problem with your input data.
            Please check your input data and try again. If you continue to have problems, please open a new issue including
            the original command and the error above:  <https://github.com/nextstrain/augur/issues/new/choose>
            """))
        sys.exit(2)
    except Exception:
        traceback.print_exc(file=sys.stderr)
        print_err("\n")
        print_err(dedent("""\
            An error occurred (see above) that has not been properly handled by Augur.
            To report this, please open a new issue including the original command and the error above:
                <https://github.com/nextstrain/augur/issues/new/choose>
            """))
        sys.exit(2)


def add_version_alias(parser):
    """
    Add --version as a (hidden) alias for the version command.

    It's not uncommon to blindly run a command with --version as the sole
    argument, so its useful to make that Just Work.
    """

    class run_version_command(argparse.Action):
        def __call__(self, *args, **kwargs):
            opts = SimpleNamespace()
            sys.exit( version.run(opts) )  # noqa: F821; This is imported from version_file.

    return parser.add_argument(
        "--version",
        nargs  = 0,
        help   = argparse.SUPPRESS,
        action = run_version_command)
