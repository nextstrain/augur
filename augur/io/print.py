import gettext
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


def indented_list(xs, prefix):
    return f"\n{prefix}".join(xs)


# Use ngettext() without a message catalog for its singular/plural handling so
# we can make proper error messages.  gettext() (no "n") is conventionally
# aliased as "_", so alias ngettext() as "_n".
_n = gettext.NullTranslations().ngettext
