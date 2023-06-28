import pytest
from pathlib import Path
from shutil import which
from subprocess import run

topdir = Path(__file__).resolve().parent.parent

pyright = which("pyright")
npx = which("npx")

if pyright:
    pyright = [pyright]
elif npx:
    pyright = [npx, "pyright"]
else:
    pyright = None

@pytest.mark.skipif(not pyright, reason = "pyright is not available")
def test_pyright():
    # Check the exit status ourselves for nicer test output on failure
    result = run(pyright, cwd = topdir)
    assert result.returncode == 0, "pyright exited with errors"
