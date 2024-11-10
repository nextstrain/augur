r"""
Write a file like Augur, with transparent optimized compression and universal newlines.

Input is always from stdin.

Output is to the given file path, as the compression format detection require
it.  The given path may be "-" to explicitly write to stdout, but no
decompression will be done.

Universal newline translation is always performed, so \n, \r\n, and \r in the
input are all translated to the system's native newlines (e.g. \n on Unix, \r\n
on Windows) in the output.
"""
import io
import sys
from shutil import copyfileobj

from .io.file import open_file
from .utils import first_line


# The buffer size used by xopen() (which underlies open_file()), which notes:
# 128KB [KiB] buffer size also used by cat, pigz etc. It is faster than the 8K
# [KiB] default.
BUFFER_SIZE = max(io.DEFAULT_BUFFER_SIZE, 128 * 1024)


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("write-file", help=first_line(__doc__))
    parser.add_argument("path", metavar="PATH", help="path to file")
    return parser


def run(args):
    with open_file(args.path, "wt", newline=None) as f:
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
        copyfileobj(sys.stdin, f, BUFFER_SIZE)
