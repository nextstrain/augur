import json

from augur.__version__ import __version__
from augur.__version__ import is_augur_version_compatible
from augur.errors import AugurError
from augur.io.file import open_file
from augur.io.print import print_err
from augur.types import ValidationMode
from augur.validate import validate_json, ValidateError, load_json_schema


FILTERED_ATTRS = ["generated_by"]


class NodeDataFile:
    def __init__(self, fname, validation_mode=ValidationMode.ERROR):
        self.fname = fname
        self.validation_mode = validation_mode

        with open_file(fname) as jfile:
            self.attrs = json.load(jfile)

        self.validate()

    @property
    def annotations(self):
        return self.attrs.get("annotations")

    @property
    def nodes(self):
        return self.attrs.get("nodes", {})

    @property
    def branches(self):
        # these are optional, so we provide an empty dict as a default
        return self.attrs.get("branches", {})

    @property
    def generated_by(self):
        return self.attrs.get("generated_by")

    @property
    def is_generated_by_incompatible_augur(self):
        if not isinstance(self.generated_by, dict):
            # If it's not a dict created by augur, we can't reliably classify it as incompatible
            return False

        generated_by_augur = self.generated_by.get("program") == "augur"
        compatible_version = is_augur_version_compatible(
            self.generated_by.get("version")
        )

        return generated_by_augur and not compatible_version

    def items(self):
        filtered_attrs = {
            key: value for key, value in self.attrs.items() if key not in FILTERED_ATTRS
        }

        return filtered_attrs.items()

    def validate(self):
        if self.validation_mode is ValidationMode.SKIP:
            return # don't perform validation (i.e. skip it)

        try:
            if self.annotations:
                try:
                    # validate_json will print any errors to STDERR
                    validate_json(
                        self.annotations,
                        load_json_schema("schema-annotations.json"),
                        self.fname,
                    )
                except ValidateError as err:
                    raise ValidateError(
                        f"{self.fname} contains an `annotations` attribute of an invalid JSON format. Was it "
                        "produced by different version of augur the one you are currently using "
                        f" ({__version__})? Please check the program that produced that JSON file."
                    ) from err

            if not isinstance(self.nodes, dict):
                raise ValidateError(
                    f"`nodes` value in {self.fname} is not a dictionary. Please check the formatting of this JSON!"
                )

            if not isinstance(self.branches, dict):
                raise ValidateError(
                    f"`branches` value in {self.fname} is not a dictionary. Please check the formatting of this JSON!"
                )

            if not self.nodes and not self.branches:
                print_err(
                    f"WARNING: {self.fname} has empty or nonexistent `nodes` and `branches`. Please check the formatting of this JSON!"
                )

            if self.is_generated_by_incompatible_augur:
                msg = (
                    f"Augur version incompatibility detected: the JSON {self.fname} was generated by "
                    f"{self.generated_by}, which is incompatible with the current augur version "
                    f"({__version__}). We suggest you rerun the pipeline using the current version of "
                    "augur."
                )
                raise ValidateError(msg)
        except ValidateError as e:
            if self.validation_mode is ValidationMode.WARN:
                # string-representation of error mirrors how augur ultimately handles AugurError
                print_err(f"WARNING: {e}")
            else:
                # Re-raising within an except block will use implicit chaining (via __context__)
                # however since we ultimately catch AugurError and use a string representation
                # this isn't particularly helpful. One day we may add a dev-flag to allow printing
                # the traceback or similar
                raise AugurError(*e.args)


class NodeDataObject(NodeDataFile):
    """
    NodeDataObject is identical to NodeDataFile except it takes a node-data dict
    rather than loading the node data from a file
    """
    def __init__(self, node_data_json, fname, validation_mode=ValidationMode.ERROR):
        self.fname = fname
        self.validation_mode = validation_mode
        self.attrs = node_data_json
        self.validate()
