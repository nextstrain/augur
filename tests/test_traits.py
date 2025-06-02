import argparse
import pytest
from augur import traits
from augur.errors import AugurError

@pytest.fixture
def branch_labels_arg():
    parser = argparse.ArgumentParser(prog="augur")
    traits.register_parser(parser.add_subparsers())
    def parse(branch_labels):
        arg_str = f"traits --tree X --metadata Y --columns A B C {branch_labels}"
        return parser.parse_args(arg_str.split(" ")).branch_labels
    return parse

class TestBranchLabellerArgParsing:
    def test_no_args(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg(""), column)
        assert bl.column == None

    def test_flag_arg(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == None

    def test_arg_for_another_column(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels B"), column)
        assert bl.column == None

    def test_arg_column_name_only(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels A"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == None

    def test_arg_column_name_then_colon(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels A:"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == None

    def test_arg_column_name_with_confidence(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels A:75"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == 75

    def test_arg_column_name_with_label_and_confidence(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels A:jumps:75"), column)
        assert bl.column == column
        assert bl.column_label == "jumps"
        assert bl._confidence_threshold == 75

    def test_arg_column_name_with_label_but_no_confidence(self, branch_labels_arg):
        column = "A"
        bl = traits.BranchLabeller(branch_labels_arg("--branch-labels A:jumps:"), column)
        assert bl.column == column
        assert bl.column_label == "jumps"
        assert bl._confidence_threshold == None

    def test_non_numeric_confidence(self, branch_labels_arg):
        column = "A"
        with pytest.raises(AugurError):
            traits.BranchLabeller(branch_labels_arg("--branch-labels A:jumps"), column)

    def test_missing_column_label(self, branch_labels_arg):
        column = "A"
        with pytest.raises(AugurError):
            traits.BranchLabeller(branch_labels_arg("--branch-labels A::75"), column)

    def test_too_many_params(self, branch_labels_arg):
        column = "A"
        with pytest.raises(AugurError):
            traits.BranchLabeller(branch_labels_arg("--branch-labels A:jumps:75:extra"), column)
