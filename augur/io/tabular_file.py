import csv
from typing import List

from augur.errors import AugurError
from .file import File


class TabularFile(File):
    """Represents a tabular file with a delimiter and optional header row.

    A pandas DataFrame can provide these properties, however this class is more
    suitable for large files as it does not store the rows in memory.

    The properties on this class can also be used on sub-classes for more
    context-specific usage.
    """

    def __init__(self, file: str, header: bool = None, names: List[str] = None):
        super().__init__(file)

        self.header = header or True
        """Indicates whether a header is present in the file."""

        self.names = names or []
        """Column names."""

        if (not self.header and not self.names) or (self.header and self.names):
            raise ValueError("Tabular file must have either a header row or column names specified, but not both.")

        self.delimiter = self._get_delimiter()
        """Delimiter of tabular file."""

        if not self.names:
            with self.open() as f:
                reader = csv.reader(f, delimiter=self.delimiter)
                try:
                    self.names = next(reader)
                except StopIteration:
                    raise AugurError(f"Expected a header row in '{self.file}' but it is empty.")

    def rows(self):
        """Yield rows from a tabular file."""
        with self.open() as f:
            reader = csv.DictReader(f, delimiter=self.delimiter, fieldnames=self.names)
            if self.header:
                next(reader)
            for row in reader:
                if not row:
                    # Skip empty lines
                    continue
                yield row

    def _get_delimiter(self):
        """Get the delimiter of a tabular file."""
        possible_delimiters = [',', '\t']

        with self.open() as f:
            try:
                # Infer delimiter from the first line.
                dialect = csv.Sniffer().sniff(f.readline(), delimiters=possible_delimiters)
            except csv.Error:
                # This can happen for single-column files, e.g. VCF sequence indexes
                # If so, use a tab character as the default.
                return '\t'
        return dialect.delimiter
