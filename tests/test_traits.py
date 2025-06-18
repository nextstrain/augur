import argparse
import pytest
from augur import traits
from augur.errors import AugurError

@pytest.fixture
def branch_labels_args():
    parser = argparse.ArgumentParser(prog="augur")
    traits.register_parser(parser.add_subparsers())
    def parse(branch_labels_args):
        arg_str = f"traits --tree X --metadata Y --columns A B C {branch_labels_args}"
        args = parser.parse_args(arg_str.split(" "))
        return [args.branch_labels, args.branch_confidence]
    return parse

class TestBranchLabellerArgParsing:
    def test_no_args(self, branch_labels_args):
        column = "A"
        bl = traits.BranchLabeller(*branch_labels_args(""), column)
        assert bl.column == None

    def test_flag_labels_arg_is_invalid(self, branch_labels_args):
        column = "A"
        with pytest.raises(BaseException):
            traits.BranchLabeller(*branch_labels_args("--branch-labels"), column)

    def test_arg_for_another_column(self, branch_labels_args):
        column = "A"
        bl = traits.BranchLabeller(*branch_labels_args("--branch-labels B"), column)
        assert bl.column == None

    def test_arg_column_name_only(self, branch_labels_args):
        column = "A"
        bl = traits.BranchLabeller(*branch_labels_args("--branch-labels A"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == None

    def test_custom_label_name(self, branch_labels_args):
        column = "A"
        bl = traits.BranchLabeller(*branch_labels_args("--branch-labels A=jumps"), column)
        assert bl.column == column
        assert bl.column_label == "jumps"
        assert bl._confidence_threshold == None

    def test_too_many_equals_labels(self, branch_labels_args):
        column = "A"
        with pytest.raises(AugurError):
            traits.BranchLabeller(*branch_labels_args("--branch-labels A=jumps=75"), column)

    def test_confidence_arg_is_invalid(self, branch_labels_args):
        column = "A"
        with pytest.raises(BaseException):
            traits.BranchLabeller(*branch_labels_args("--branch-confidence"), column)

    def test_confidence_arg(self, branch_labels_args):
        column = "A"
        bl = traits.BranchLabeller(*branch_labels_args("--branch-labels A --branch-confidence A=75"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == 75

    def test_confidence_arg_for_different_column(self, branch_labels_args):
        column = "A"
        bl = traits.BranchLabeller(*branch_labels_args("--branch-labels A --branch-confidence B=75"), column)
        assert bl.column == column
        assert bl.column_label == column
        assert bl._confidence_threshold == None

    def test_non_numeric_confidence(self, branch_labels_args):
        column = "A"
        with pytest.raises(AugurError):
            traits.BranchLabeller(*branch_labels_args("--branch-labels A --branch-confidence A=foo"), column)

    def test_too_many_equals_confidence(self, branch_labels_args):
        column = "A"
        with pytest.raises(AugurError):
            traits.BranchLabeller(*branch_labels_args("--branch-labels A --branch-confidence A=B=75"), column)
