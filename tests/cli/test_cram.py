import pytest
import pathlib
import cram


# Locate all Cram tests to be executed.
test_dir = "./tests/cli"
cram_tests = list(pathlib.Path(test_dir).glob("**/*.t"))


def get_ids(args):
    """
    Convert PosixPaths representing individual Cram test files to string.
    This allows Pytest to include each Cram test's filepath and success/failure status in the pytest output.
    """
    return str(args)


class TestCram:
    """
    Pytest wrapper for Cram CLI tests.
    This test executes each Cram test defined under tests/cli/, and fails if Cram returns any diffs.

    To generate a coverage report for all Cram CLI tests, execute this command:
        pytest -vv -s tests/cli/test_cram.py --cov augur --cov-report term
    """

    @pytest.mark.xfail(
        reason="Cram tests currently fail for export, filter, mask, refine, and tree."
    )
    @pytest.mark.parametrize("cram_test_file", cram_tests, ids=get_ids)
    def test_all(self, cram_test_file):
        # cram expects a bytes literal here, e.g. b"tests/cli/ancestral/add_to_alignment/ancestral/ancestral.t"
        ins, outs, diffs = cram.testfile(path=bytes(cram_test_file))
        diff_list = list(diffs)
        assert len(diff_list) == 0, f"Cram diffs: {diff_list}"
