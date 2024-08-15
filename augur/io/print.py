import sys

from augur.debug import DEBUGGING


def print_err(*args):
    """Print to stderr. When data goes to stdout (most cases), this should be
    used for any informational messages, not just errors/warnings."""
    print(*args, file=sys.stderr)


def print_debug(*args):
    """Print to stderr if in debugging mode."""
    if DEBUGGING:
        print_err(*args)
