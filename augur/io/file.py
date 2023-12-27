import os
from contextlib import contextmanager
from io import IOBase
from textwrap import dedent
from xopen import PipedCompressionReader, PipedCompressionWriter, xopen
from augur.errors import AugurError


ENCODING = "utf-8"


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
            raise AugurError(dedent(f"""\
                File '{path_or_buffer}' contains {e.object[e.start:e.end]} which is not encoded as '{e.encoding}'.
                Try saving the file using '{e.encoding}' encoding."""))


    elif isinstance(path_or_buffer, (IOBase, PipedCompressionReader, PipedCompressionWriter)):
        yield path_or_buffer

    else:
        raise TypeError(f"Type {type(path_or_buffer)} is not supported.")
