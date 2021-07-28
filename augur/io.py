#!/usr/bin/env python3
"""Interfaces for reading and writing data also known as input/output (I/O)
"""
import Bio.SeqIO
import Bio.SeqRecord
from contextlib import contextmanager
import pandas as pd
from pathlib import Path
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


def read_metadata(metadata_file, valid_index_cols=("strain", "name"), **kwargs):
    """Read metadata from a given filename and into a pandas `DataFrame` or
    `TextFileReader` object.

    Parameters
    ----------
    metadata_file : str
        Path to a metadata file to load.
    valid_index_cols : list[str]
        List of possible index column names to check for, ordered by priority.
    kwargs : dict
        Keyword arguments to pass through to pandas.read_csv.

    Returns
    -------
    pandas.DataFrame or pandas.TextFileReader

    Raises
    ------
    KeyError :
        When the metadata file does not have any valid index columns.

    """
    default_kwargs = {
        "sep": None,
        "engine": "python",
        "skipinitialspace": True,
        "dtype": {
            "strain": "string",
            "name": "string,"
        }
    }
    default_kwargs.update(kwargs)

    # Inspect the first chunk of the metadata, to find any valid index columns.
    chunk = pd.read_csv(
        metadata_file,
        iterator=True,
        **default_kwargs,
    ).read(nrows=1)

    index_cols = [
        valid_index_col
        for valid_index_col in valid_index_cols
        if valid_index_col in chunk.columns
    ]

    # If we couldn't find a valid index column in the metadata, alert the user.
    if len(index_cols) == 0:
        raise KeyError(f"Could not find any valid index columns from the list of possible columns ({valid_index_cols})")
    else:
        index_col = index_cols[0]

    default_kwargs["index_col"] = index_col
    return pd.read_csv(
        metadata_file,
        **default_kwargs
    )


def read_sequences(*paths, format="fasta"):
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc.) and return a stream
    of sequence records in the requested format (e.g., "fasta", "genbank", etc.).

    Parameters
    ----------
    paths : list of str or Path-like objects
        One or more paths to sequence files of any type supported by BioPython.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.).

    Yields
    ------
    Bio.SeqRecord.SeqRecord
        Sequence record from the given path(s).

    """
    for path in paths:
        # Open the given path as a handle, inferring the file's compression.
        # This way we can pass a handle to BioPython's SeqIO interface
        # regardless of the compression mode.
        with open_file(path) as handle:
            sequences = Bio.SeqIO.parse(handle, format)

            for sequence in sequences:
                yield sequence


def write_sequences(sequences, path_or_buffer, format="fasta"):
    """Write sequences to a given path in the given format.

    Automatically infer compression mode (e.g., gzip, etc.) based on the path's
    filename extension.

    Parameters
    ----------
    sequences : iterable of Bio.SeqRecord.SeqRecord objects
        A list-like collection of sequences to write

    path_or_buffer : str or Path-like object or IO buffer
        A path to a file to write the given sequences in the given format.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.)

    Returns
    -------
    int :
        Number of sequences written out to the given path.

    """
    with open_file(path_or_buffer, "wt") as handle:
        # Bio.SeqIO supports writing to the same handle multiple times for specific
        # file formats. For the formats we use, this function call should work for
        # both a newly opened file handle or one that is provided by the caller.
        # For more details see:
        # https://github.com/biopython/biopython/blob/25f5152f4aeefe184a323db25694fbfe0593f0e2/Bio/SeqIO/__init__.py#L233-L251
        sequences_written = Bio.SeqIO.write(
            sequences,
            handle,
            format
        )

    return sequences_written
