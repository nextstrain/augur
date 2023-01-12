from contextlib import contextmanager
from xopen import xopen


@contextmanager
def open_file(path_or_buffer, mode="r", **kwargs):
    """Opens a given file path and returns the handle.

    Transparently handles compressed inputs and outputs.

    Parameters
    ----------
    path_or_buffer : str or Path-like or IO buffer
        Name of the file to open or an existing IO buffer

    mode : str
        Mode to open file (read or write)

    Returns
    -------
    IO
        File handle object

    """
    try:
        with xopen(path_or_buffer, mode, **kwargs) as handle:
            yield handle
    except TypeError:
        yield path_or_buffer


class File:
    """Represents a file."""

    def __init__(self, file: str):
        self.file = file
        """File to be loaded."""

    def open(self, **kwargs):
        """Open the file with auto-compression/decompression."""
        return open_file(self.file, **kwargs)
