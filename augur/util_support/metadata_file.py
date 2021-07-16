import functools
import pandas
import sys


class MetadataFile:
    """
    Represents a CSV or TSV file containing metadata

    The file must contain exactly one of a column named `strain` or `name`,
    which is used to match metadata with samples.
    """

    def __init__(self, fname, query=None, as_data_frame=False):
        self.fname = fname
        self.query = query
        self.as_data_frame = as_data_frame

        self.key_type = self.find_key_type()

    def read(self):
        self.check_metadata_duplicates()

        # augur assumes the metadata dict will contain either "strain" or "name" (the
        # indexed column), but DataFrame.to_dict("index") does not place the indexed
        # column in the dict. So let's make a copy of the indexed column so that the
        # original "strain"/"name" remains in the output.
        self.metadata["_index"] = self.metadata[self.key_type]

        metadata = self.metadata.set_index("_index")

        if self.as_data_frame:
            return metadata, self.columns
        else:
            return metadata.to_dict("index"), self.columns

    @property
    @functools.lru_cache()
    def metadata(self):
        """
        Return list of dicts representing the metadata in the file.

        If a query was supplied, apply it.
        """

        metadata = self.parse_file()

        if self.query:
            try:
                metadata = metadata.query(self.query).copy()
            except Exception as e:
                raise ValueError(
                    f"Error applying pandas query to metadata: `{self.query}` ({e})"
                )

        return metadata

    def check_metadata_duplicates(self):
        duplicates = (
            self.metadata[self.key_type]
            .value_counts()
            .reset_index()
            .query(f"{self.key_type} > 1")["index"]
            .values
        )

        if len(duplicates) > 0:
            raise ValueError(
                f"Duplicated {self.key_type} in metadata: {', '.join(duplicates)}"
            )

    @property
    @functools.lru_cache()
    def columns(self):
        return list(self.parse_file().columns)

    def find_key_type(self):
        if "strain" not in self.columns and "name" not in self.columns:
            raise ValueError(
                f"Metadata file {self.fname} does not contain `name` or `strain`"
            )

        if "strain" in self.columns and "name" in self.columns:
            print(
                f"WARNING: Metadata file {self.fname} contains both `name` and `strain`. Using `strain`.",
                file=sys.stderr,
            )

        if "strain" in self.columns:
            return "strain"

        return "name"

    @functools.lru_cache()
    def parse_file(self):
        return pandas.read_csv(
            self.fname,
            sep=None,  # csv.Sniffer will automatically detect sep
            engine="python",
            skipinitialspace=True,
            dtype={"strain":"string", "name":"string"}
        ).fillna("")
