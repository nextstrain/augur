import pytest
import sys
import os
from shutil import copyfile

# Following https://stackoverflow.com/questions/13493288/python-cli-program-unit-testing/13500346#13500346
# and https://docs.pytest.org/en/latest/reference.html#testdir

pytest_plugins = "pytester"

@pytest.fixture
def run(testdir):
    def do_run(*args):
        args = ["augur", "tree"] + list(args)
        return testdir.run(*args)
    return do_run

def is_valid_newick(contents):
    return True

def test_looks_like_a_tree(tmpdir, run):
    copyfile(os.path.join(sys.path[0], "./assets/aligned.fasta"), tmpdir.join("aligned.fasta"))
    output = tmpdir.join("tree.new")
    result = run("--alignment", tmpdir.join("aligned.fasta"), "--method", "iqtree", "--output", output)
    assert result.ret == 0
    with output.open("r") as f:
        newick = f.read()
    assert(is_valid_newick(newick))
