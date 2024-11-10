import os
from contextlib import contextmanager
from io import IOBase
from textwrap import dedent
from xopen import xopen
from augur.errors import AugurError

# Workaround to maintain compatibility with both xopen v1 and v2
# Around November 2024, we shall drop support for xopen v1
# by removing the try-except block and using
# _PipedCompressionProgram directly
try:
    from xopen import _PipedCompressionProgram as PipedCompressionReader
    from xopen import _PipedCompressionProgram as PipedCompressionWriter
except ImportError:
    from xopen import (  # type: ignore[attr-defined, no-redef]  
        PipedCompressionReader,
        PipedCompressionWriter,
    )

ENCODING = "utf-8"

PANDAS_READ_CSV_OPTIONS = {
    'encoding': ENCODING,
}


@contextmanager
def open_file(path_or_buffer, mode="r", **kwargs):
    """Opens a given file path and returns the handle.

    Transparently handles compressed inputs and outputs.

    Parameters
    ----------
    path_or_buffer
        Name of the file to open or an existing IO buffer

    mode : str
        Mode to open file (read or write)

    Returns
    -------
    IO
        File handle object

    """

    # Read all files using a specific encoding.
    kwargs['encoding'] = ENCODING

    if isinstance(path_or_buffer, (str, os.PathLike)):
        try:
            with xopen(path_or_buffer, mode, **kwargs) as handle:
                yield handle
        except UnicodeDecodeError as e:
            # TODO: Consider moving this to the top-level error handler to
            # handle errors from other I/O functions such as pandas.read_csv.
            # This is not trivial since the filepath is useful to include in the
            # message, but is not available through UnicodeDecodeError alone.
            raise AugurError(dedent(f"""\
                File {path_or_buffer!r} contains {e.object[e.start:e.end]!r} which is not valid in the expected {e.encoding!r} encoding.
                Try re-saving the file using the {e.encoding!r} encoding."""))


    elif isinstance(path_or_buffer, (IOBase, PipedCompressionReader, PipedCompressionWriter)):
        yield path_or_buffer

    else:
        raise TypeError(f"Type {type(path_or_buffer)} is not supported.")
