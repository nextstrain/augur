import Bio.Phylo
from io import StringIO
from pathlib import Path
import pytest
import sys

# we assume (and assert) that this script is running from the tests/ directory
sys.path.append(str(Path(__file__).parent.parent.parent))

from augur.export_v2 import convert_tree_to_json_structure
from augur.validate import ValidateError
from augur.validate_export import ensure_no_duplicate_names


class TestValidateExport():
    def test_export_without_duplicate_names(self):
        # Create a tree with unique tip names.
        tree = Bio.Phylo.read(StringIO("root(A, internal(B, C))"), "newick")
        metadata = {"A": {}, "B": {}, "C": {}, "root": {}, "internal": {}}
        root = convert_tree_to_json_structure(tree.root, metadata, None)
        ensure_no_duplicate_names(root, ValidateError)

    def test_export_with_duplicate_names(self):
        # Create a tree with duplicate tip names.
        tree = Bio.Phylo.read(StringIO("root(A, internal(B, B))"), "newick")
        metadata = {"A": {}, "B": {}, "root": {}, "internal": {}}
        root = convert_tree_to_json_structure(tree.root, metadata, None)

        with pytest.raises(ValidateError):
            ensure_no_duplicate_names(root, ValidateError)
