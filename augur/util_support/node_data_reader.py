import Bio.Phylo
import sys

from augur.types import ValidationMode
from augur.util_support.node_data import DuplicatedNonDictAttributeError
from augur.util_support.node_data import NodeData
from augur.util_support.node_data_file import NodeDataFile


class NodeDataReader:
    """
    Parses one or more node data files and combines them into a single data structure.

    The format of Node data files is described in `docs/usage/json_format.md`. The `generated_by` attribute is
    not included in the resulting data structure.

    If a tree file is specified, it is used to verify the node names.

    If validation_mode is set to :py:attr:`augur.types.ValidationMode.SKIP` no validation is performed.
    """

    def __init__(self, filenames, tree_file=None, validation_mode=ValidationMode.ERROR):
        if not isinstance(filenames, list):
            filenames = [filenames]
        self.filenames = filenames
        self.tree_file = tree_file
        self.validation_mode = validation_mode

    def read(self):
        node_data = self.build_node_data()

        self.check_against_tree_file(node_data)

        return node_data

    def build_node_data(self):
        node_data = NodeData()

        for node_data_file in self.node_data_files:
            try:
                node_data.update(node_data_file)
            except DuplicatedNonDictAttributeError as e:
                print(
                    f"Error parsing {node_data_file.fname}: Mixed dict and non-dict values found for "
                    f"attribute `{e.key}`. Augur requires attributes of the same name to be uniformly dicts "
                    "or non-dicts.",
                    file=sys.stderr,
                )
                sys.exit(2)

        return node_data.attrs

    @property
    def node_data_files(self):
        return (NodeDataFile(fname, validation_mode = self.validation_mode) for fname in self.filenames)

    def check_against_tree_file(self, node_data):
        if not self.tree_file:
            return

        if set(node_data["nodes"].keys()) != self.node_names_from_tree_file:
            print(
                f"Names of nodes (including internal nodes) of tree {self.tree_file} don't match node names "
                "in the node data files.",
                file=sys.stderr,
            )
            sys.exit(2)

    @property
    def node_names_from_tree_file(self):
        try:
            tree = Bio.Phylo.read(self.tree_file, "newick")
        except Exception as e:
            print(
                f"Failed to read tree from file {self.tree_file}: {e}", file=sys.stderr
            )
            sys.exit(2)

        return set([clade.name for clade in tree.find_clades()])
