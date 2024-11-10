r"""
Read a file like Augur, with transparent optimized decompression and universal newlines.

Input is read from the given file path, as the compression format detection
requires a seekable stream.  The given path may be "-" to explicitly read from
stdin, but no decompression will be done.

Output is always to stdout.

Universal newline translation is always performed, so \n, \r\n, and \r in the
input are all translated to the system's native newlines (e.g. \n on Unix, \r\n
on Windows) in the output.
"""
import io
import os
import signal
import sys
from shutil import copyfileobj

from .io.file import open_file
from .utils import first_line


# The buffer size used by xopen() (which underlies open_file()), which notes:
# 128KB [KiB] buffer size also used by cat, pigz etc. It is faster than the 8K
# [KiB] default.
BUFFER_SIZE = max(io.DEFAULT_BUFFER_SIZE, 128 * 1024)

SIGPIPE = getattr(signal, "SIGPIPE", None)


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("read-file", help=first_line(__doc__))
    parser.add_argument("path", metavar="PATH", help="path to file")
    return parser


def run(args):
    with open_file(args.path, "rt", newline=None) as f:
        # It's tempting to want to splice(2) here, but it turns out to make
        # little sense.  Firstly, the availability of splice(2) is Linux,
        # Python ≥3.10, and one of the files needs to be a pipe.  The chance of
        # all of those together is slim-to-none, particularly because even in
        # the common case of xopen() reading from a pipe—the stdout of an
        # external decompression process—we can't use that pipe directly
        # because xopen() always buffers the first block of the file into
        # Python.¹  Secondly, we want universal newline handling—so that
        # callers get behaviour consistent with the rest of Augur—and that
        # rules out splice(2).
        #
        # Copying the data thru Python instead of with splice(2) seems fast
        # enough in some quick trials with large files (e.g. against `zstd`
        # directly), and the bottlenecks in pipelines using this command will
        # often not be this command's i/o.
        #   -trs, 11 July 2024, updated 24 July 2024
        #
        # ¹ <https://github.com/pycompression/xopen/blob/67651844/src/xopen/__init__.py#L374>

        # Handle SIGPIPE, which Python converts to BrokenPipeError, gracefully
        # and like most Unix programs.  See also
        # <https://docs.python.org/3/library/signal.html#note-on-sigpipe>.
        try:
            copyfileobj(f, sys.stdout, BUFFER_SIZE)

            # Force a flush so if SIGPIPE is going to happen it happens now.
            sys.stdout.flush()
        except BrokenPipeError:
            # Avoid errors from Python automatically flushing stdout on exit.
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())

            # Return conventional exit status for "killed by SIGPIPE" on Unix.
            return 128 + SIGPIPE if SIGPIPE else 1
