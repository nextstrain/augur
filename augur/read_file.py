r"""
Read one or more files like Augur, with transparent optimized decompression and
universal newlines. Supported compression formats: gzip (.gz), bzip2 (.bz2),
xz (.xz), zstandard (.zst).

Input is read from each given file path, as the compression format detection
requires a seekable stream. A path may be "-" to explicitly read from stdin, but
no decompression will be done.

Output from each file is concatenated together and written to stdout.

Universal newline translation is always performed, so \n, \r\n, and \r in the
input are all translated to the system's native newlines (e.g. \n on Unix, \r\n
on Windows) in the output. Additionally, each file is standardized to have
trailing newlines.
"""
import io
import os
import signal
import sys

from .argparse_ import ExtendOverwriteDefault
from .io.file import open_file
from .utils import first_line


# The buffer size used by xopen() (which underlies open_file()), which notes:
# 128KB [KiB] buffer size also used by cat, pigz etc. It is faster than the 8K
# [KiB] default.
BUFFER_SIZE = max(io.DEFAULT_BUFFER_SIZE, 128 * 1024)

SIGPIPE = getattr(signal, "SIGPIPE", None)


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("read-file", help=first_line(__doc__))
    parser.add_argument("paths", metavar="PATH", nargs="+", action=ExtendOverwriteDefault, help="paths to files")
    return parser


def run(args):
    for path in args.paths:
        with open_file(path, "rt", newline=None) as f:
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
                copyfileobj(f, sys.stdout, ensure_trailing_newline=True)

                # Force a flush so if SIGPIPE is going to happen it happens now.
                sys.stdout.flush()
            except BrokenPipeError:
                # Avoid errors from Python automatically flushing stdout on exit.
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())

                # Return conventional exit status for "killed by SIGPIPE" on Unix.
                return 128 + SIGPIPE if SIGPIPE else 1


def copyfileobj(fsrc, fdst, length=0, ensure_trailing_newline=False):
    """copy data from file-like object fsrc to file-like object fdst"""
    if not length:
        length = BUFFER_SIZE

    # Localize variable access to minimize overhead.
    fsrc_read = fsrc.read
    fdst_write = fdst.write

    last_char = None
    while buf := fsrc_read(length):
        fdst_write(buf)
        last_char = buf[-1]

    if ensure_trailing_newline and last_char != '\n':
        fdst_write('\n')
