import os
from contextlib import contextmanager
from io import IOBase
from textwrap import dedent
from xopen import xopen, _PipedCompressionProgram
from augur.errors import AugurError


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


    elif isinstance(path_or_buffer, (IOBase, _PipedCompressionProgram)):
        yield path_or_buffer

    else:
        raise TypeError(f"Type {type(path_or_buffer)} is not supported.")
