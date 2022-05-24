import pandas as pd


def read_metadata(metadata_file, id_columns=("strain", "name"), chunk_size=None):
    """Read metadata from a given filename and into a pandas `DataFrame` or
    `TextFileReader` object.

    Parameters
    ----------
    metadata_file : str
        Path to a metadata file to load.
    id_columns : list[str]
        List of possible id column names to check for, ordered by priority.
    chunk_size : int
        Size of chunks to stream from disk with an iterator instead of loading the entire input file into memory.

    Returns
    -------
    pandas.DataFrame or pandas.TextFileReader

    Raises
    ------
    KeyError :
        When the metadata file does not have any valid index columns.


    For standard use, request a metadata file and get a pandas DataFrame.

    >>> read_metadata("tests/functional/filter/data/metadata.tsv").index.values[0]
    'COL/FLR_00024/2015'

    Requesting an index column that doesn't exist should produce an error.

    >>> read_metadata("tests/functional/filter/data/metadata.tsv", id_columns=("Virus name",))
    Traceback (most recent call last):
      ...
    Exception: None of the possible id columns (('Virus name',)) were found in the metadata's columns ('strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'city', 'db', 'segment', 'authors', 'url', 'title', 'journal', 'paper_url')

    We also allow iterating through metadata in fixed chunk sizes.

    >>> for chunk in read_metadata("tests/functional/filter/data/metadata.tsv", chunk_size=5):
    ...     print(chunk.shape)
    ...
    (5, 14)
    (5, 14)
    (2, 14)

    """
    kwargs = {
        "sep": None,
        "engine": "python",
        "skipinitialspace": True,
        "na_filter": False,
    }

    if chunk_size:
        kwargs["chunksize"] = chunk_size

    # Inspect the first chunk of the metadata, to find any valid index columns.
    metadata = pd.read_csv(
        metadata_file,
        iterator=True,
        **kwargs,
    )
    chunk = metadata.read(nrows=1)
    metadata.close()

    id_columns_present = [
        id_column
        for id_column in id_columns
        if id_column in chunk.columns
    ]

    # If we couldn't find a valid index column in the metadata, alert the user.
    if not id_columns_present:
        raise Exception(f"None of the possible id columns ({id_columns!r}) were found in the metadata's columns {tuple(chunk.columns)!r}")
    else:
        index_col = id_columns_present[0]

    # If we found a valid column to index the DataFrame, specify that column and
    # also tell pandas that the column should be treated like a string instead
    # of having its type inferred. This latter argument allows users to provide
    # numerical ids that don't get converted to numbers by pandas.
    kwargs["index_col"] = index_col
    kwargs["dtype"] = {index_col: "string"}

    return pd.read_csv(
        metadata_file,
        **kwargs
    )
