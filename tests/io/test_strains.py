from pathlib import Path

from augur.io.strains import read_strains


def test_read_strains(tmpdir):
    # Write one list of filenames with some unnecessary whitespace.
    strains1 = Path(tmpdir) / Path("strains1.txt")
    with open(strains1, "w") as oh:
        oh.write("strain1 # this is an inline comment about strain 1\nstrain2\n   # this is a comment preceded by whitespace.\n")

    # Write another list of filenames with a comment.
    strains2 = Path(tmpdir) / Path("strains2.txt")
    with open(strains2, "w") as oh:
        oh.write("# this is a comment. ignore this.\nstrain2\nstrain3\n")

    strains = read_strains(strains1, strains2)
    assert len(strains) == 3
    assert "strain1" in strains
