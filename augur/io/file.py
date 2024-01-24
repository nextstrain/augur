import os
from contextlib import contextmanager
from io import IOBase
from xopen import PipedCompressionReader, PipedCompressionWriter, xopen


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
        with xopen(path_or_buffer, mode, **kwargs) as handle:
            yield handle

    elif isinstance(path_or_buffer, (IOBase, PipedCompressionReader, PipedCompressionWriter)):
        yield path_or_buffer

    else:
        raise TypeError(f"Type {type(path_or_buffer)} is not supported.")
