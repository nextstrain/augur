from pathlib import Path
from subprocess import run

topdir = Path(__file__).resolve().parent.parent

def test_flake8():
    # Check the exit status ourselves for nicer test output on failure
    result = run(["flake8"], cwd=topdir, capture_output=True, text=True)
    output = result.stderr + result.stdout
    assert result.returncode == 0, f"flake8 exited with errors:\n{output}"
