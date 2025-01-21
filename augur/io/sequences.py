import Bio.SeqIO
import os

from augur.errors import AugurError
from importlib.metadata import version as installed_version
from packaging.version import Version
from typing import Iterator, Iterable, Union
from .file import open_file


def get_biopython_format(augur_format: str) -> str:
    """Validate sequence file format and return the inferred Biopython format."""
    supported_formats = {"fasta", "genbank"}
    if augur_format not in supported_formats:
        raise AugurError(f"Sequence file format {augur_format!r} is invalid. Must be one of {', '.join(repr(f) for f in sorted(supported_formats))}.")

    # Supported formats are valid Biopython formats
    biopython_format = augur_format

    # Allow comments in FASTA format using fasta-pearson in later biopython versions
    if Version(installed_version("biopython")) >= Version("1.85") and biopython_format == "fasta":
        biopython_format = "fasta-pearson"

    return biopython_format


def read_sequences(
    *paths: Iterable[Union[str, os.PathLike]],
    format: str = "fasta",
) -> Iterator[Bio.SeqIO.SeqRecord]:
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc.) and return a stream
    of sequence records given the file format.

    Parameters
    ----------
    paths
        One or more paths to sequence files.

    format
        Format of input sequences. Either "fasta" or "genbank".

    Returns
    -------
        Sequence records from the given path(s).
    """
    for path in paths:
        # Open the given path as a handle, inferring the file's compression.
        # This way we can pass a handle to BioPython's SeqIO interface
        # regardless of the compression mode.
        with open_file(path) as handle:
            sequences = Bio.SeqIO.parse(handle, get_biopython_format(format))

            for sequence in sequences:
                yield sequence


def read_single_sequence(
    path: Union[str, os.PathLike],
    format: str = "fasta",
) -> Bio.SeqIO.SeqRecord:
    """Read a single sequence from a path.

    Automatically infers compression mode.

    Parameters
    ----------
    path
        Path to a sequence file.

    format
        Format of input file. Either "fasta" or "genbank".

    Returns
    -------
        A single sequence record from the given path.
    """
    with open_file(path) as handle:
        return Bio.SeqIO.read(handle, get_biopython_format(format))


def write_sequences(sequences, path_or_buffer, format="fasta"):
    """Write sequences to a given path in the given format.

    Automatically infer compression mode (e.g., gzip, etc.) based on the path's
    filename extension.

    Parameters
    ----------
    sequences : iterable of Bio.SeqRecord.SeqRecord
        A list-like collection of sequences to write

    path_or_buffer : str or `os.PathLike` or `io.StringIO`
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


def write_records_to_fasta(records, fasta, seq_id_field='strain', seq_field='sequence'):
    """
    Write sequences from dict *records* to a *fasta* file.
    Yields the records with the *seq_field* dropped so that they can be consumed downstream.

    Parameters
    ----------
    records: iterable of dict
        Iterator that yields dict that contains sequences

    fasta: str
        Path to FASTA file

    seq_id_field: str, optional
        Field name for the sequence identifier

    seq_field: str, optional
        Field name for the genomic sequence

    Yields
    ------
    dict:
        A copy of the record with *seq_field* dropped

    Raises
    ------
    AugurError
        When the sequence id field or sequence field does not exist in a record
    """
    with open_file(fasta, "w") as output_fasta:
        for record in records:
            if seq_id_field not in record:
                raise AugurError(f"Provided sequence identifier field {seq_id_field!r} does not exist.")
            if seq_field not in record:
                raise AugurError(f"Provided sequence field {seq_field!r} does not exist.")

            output_fasta.writelines([
                f">{record[seq_id_field]}\n",
                f"{record[seq_field]}\n"
            ])

            yield {
                key: value
                for key, value in record.items()
                if key != seq_field
            }
