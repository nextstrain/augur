import csv
import pandas as pd
from typing import List

from augur.io import open_file


def get_metadata_id_column(metadata_file:str, id_columns:List[str]):
    """Returns the first column in `id_columns` that is present in the metadata.

    Raises a `ValueError` when none of `id_columns` are found.
    """
    metadata_columns = _get_column_names(metadata_file)
    for col in id_columns:
        if col in metadata_columns:
            return col
    raise ValueError(f"None of the possible id columns ({id_columns!r}) were found in the metadata's columns {tuple(metadata_columns)!r}")


def _get_column_names(file:str):
    """Get column names using pandas."""
    delimiter = get_delimiter(file)
    read_csv_kwargs = {
        "sep": delimiter,
        "engine": "c",
        "skipinitialspace": True,
        "dtype": 'string',
    }
    row = pd.read_csv(
        file,
        nrows=1,
        **read_csv_kwargs,
    )
    return list(row.columns)


def get_delimiter(file:str, delimiters:List[str]=[',', '\t']):
    """Infer tabular delimiter from first line of a file."""
    with open_file(file) as f:
        try:
            dialect = csv.Sniffer().sniff(f.readline(), delimiters=delimiters)
        except csv.Error:
            # This can happen for single-column files, e.g. VCF sequence indexes
            # If so, use a tab character as the default.
            return '\t'
    return dialect.delimiter


class TabularFileLoaderBase:
    """Base class for loading tabular files into a database.

    Stores commonly used attributes - see __init__.
    """
    def __init__(self, file:str, header=True, names:List[str]=[]):
        self.file = file
        "Path to the tabular file to be represented by an instance of this class."
        self.header = header
        "True if file has a header, otherwise false."
        self.names = names
        "List of column names in file."

        self.delimiter = get_delimiter(self.file)
        "Inferred delimiter of file."
        self.df_head = self._get_pd_df_head(nrows=100)
        "pandas.DataFrame of first 100 rows in file, used for schema purposes."

    def _get_pd_df_head(self, nrows:int) -> pd.DataFrame:
        """Get a pandas DataFrame representation of the first `nrows` rows in the file.
        
        This is used later to:
        1. Define column names
        2. Infer data types for table creation
        """
        read_csv_kwargs = {
            "sep": self.delimiter,
            "engine": "c",
            "skipinitialspace": True,
        }
        if (not self.header and not self.names) or (self.header and self.names):
            raise ValueError("DataFrame must have either a header row or column names specified, but not both.")
        if self.names:
            read_csv_kwargs['header'] = None
            read_csv_kwargs['names'] = self.names
        return pd.read_csv(
            self.file,
            nrows=nrows,
            **read_csv_kwargs,
        )

    def _iter_indexed_rows(self):
        """Yield rows from a tabular file with an additional first column for row number."""
        with open_file(self.file) as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            if self.header:
                next(reader)
            for i, row in enumerate(reader):
                if not row:
                    continue
                yield [i] + row
