import csv
from typing import Iterable, Sequence

from augur.errors import AugurError
from .file import File, open_file


DEFAULT_DELIMITERS = (',', '\t')


class InvalidDelimiter(Exception):
    pass


class TabularFile(File):
    """Represents a tabular file with a delimiter and optional header line.

    A pandas DataFrame can provide these properties, however this class is more
    suitable for large files as it does not store the rows in memory.

    The properties on this class can also be used on sub-classes for more
    context-specific usage.
    """

    columns: Sequence[str]
    """Columns extracted from the first row in the file."""

    delimiter: str
    """Inferred delimiter of file."""

    header: bool
    """Whether the first of the file represents column names."""

    def __init__(self, path: str, delimiter: str = None, delimiters: Sequence[str] = None,
                 header: bool = True, columns: Sequence[str] = None):
        """
        Parameters
        ----------
        path
            Path of the tabular file.
        delimiter
            Use this as the delimiter.
        delimiters
            Possible delimiters to use, in order of precedence.
        header
            If true, the first line will be used as column names and not a row of data.
        columns
            If set, this will be used for column names.
        """
        super().__init__(path)

        if delimiter:
            # Delimiter is given as an argument.
            if delimiters:
                raise ValueError("At most one of delimiter and delimiters can be set.")
            self.delimiter = delimiter
        else:
            # Delimiter is inferred from possible delimiters.
            if not delimiters:
                delimiters = DEFAULT_DELIMITERS

            # Infer the dialect.
            self.delimiter = get_delimiter(self.path, delimiters)

        if header:
            # Infer column names from the header.
            if columns:
                raise ValueError("Tabular file must have either a header row or column names specified, but not both.")
            else:
                with self.open() as f:
                    reader = csv.reader(f, delimiter=self.delimiter)
                    try:
                        self.header = True
                        self.columns = next(reader)
                    except StopIteration:
                        raise AugurError(f"{self.path}: Expected a header row but it is empty.")
        else:
            # Column names should be given as an argument.
            if not columns:
                raise ValueError("Tabular file must have either a header row or column names specified.")
            else:
                self.header = False
                self.columns = columns

    def rows(self, strict: bool = True):
        """Yield rows in a dictionary format. Empty lines are ignored.

        Parameters
        ----------
        strict
            If True, raise an error when a row contains more or less than the number of expected columns.
        """
        with self.open() as f:
            reader = csv.DictReader(f, delimiter=self.delimiter, fieldnames=self.columns, restkey=None, restval=None)

            if self.header:
                next(reader)

            # NOTE: Empty lines are ignored by csv.DictReader.
            # <https://github.com/python/cpython/blob/647b6cc7f16c03535cede7e1748a58ab884135b2/Lib/csv.py#L181-L185>
            for row in reader:
                if strict:
                    if None in row.keys():
                        raise AugurError(f"{self.path}: Line {reader.line_num} contains at least one extra column. The inferred delimiter is {self.delimiter!r}.")
                    if None in row.values():
                        # This is distinct from a blank value (empty string).
                        raise AugurError(f"{self.path}: Line {reader.line_num} is missing at least one column. The inferred delimiter is {self.delimiter!r}.")
                yield row


def get_delimiter(path: str, valid_delimiters: Iterable[str]):
    """Get the delimiter of a file given a list of valid delimiters."""

    for delimiter in valid_delimiters:
        if len(delimiter) != 1:
            raise AugurError(f"Delimiters must be single-character strings. {delimiter!r} does not satisfy that condition.")

    with open_file(path) as file:
        try:
            # Infer the delimiter from the first line.
            return csv.Sniffer().sniff(file.readline(), "".join(valid_delimiters)).delimiter
        except csv.Error as error:
            # This assumes all csv.Errors imply a delimiter issue. That might
            # change in a future Python version.
            raise InvalidDelimiter from error
