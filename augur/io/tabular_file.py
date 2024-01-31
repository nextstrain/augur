import csv
from typing import Iterable, Sequence

from augur.errors import AugurError
from .file import File, open_file


DEFAULT_DELIMITERS = (',', '\t')


class InvalidDelimiter(Exception):
    pass


class TabularFile(File):
    """Represents a tabular file with a delimiter and required header line.

    A pandas DataFrame can provide these properties, however this class is more
    suitable for large files as it does not store the rows in memory.

    The properties on this class can also be used on sub-classes for more
    context-specific usage.
    """

    columns: Sequence[str]
    """Columns extracted from the first row in the file."""

    delimiter: str
    """Inferred delimiter of file."""

    def __init__(self, path: str, delimiters: Sequence[str] = None):
        """
        Parameters
        ----------
        path
            Path of the tabular file.
        delimiters
            Possible delimiters to use, in order of precedence.
        """
        super().__init__(path)

        if delimiters is None:
            delimiters = DEFAULT_DELIMITERS

        # Infer the dialect.
        self.delimiter = get_delimiter(self.path, delimiters)

        # Infer the column names.
        with self.open() as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            try:
                self.columns = next(reader)
            except StopIteration:
                raise AugurError(f"{self.path}: Expected a header row but it is empty.")


    def rows(self, strict: bool = True):
        """Yield rows in a dictionary format. Empty lines are ignored.

        Parameters
        ----------
        strict
            If True, raise an error when a row contains more or less than the number of expected columns.
        """
        with self.open() as f:
            reader = csv.DictReader(f, delimiter=self.delimiter, fieldnames=self.columns, restkey=None, restval=None)

            # Skip the header row.
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
