"""
Unit tests for building a tree
"""
import sys
from pathlib import Path
import pytest
import Bio
import pandas as pd

# we assume (and assert) that this script is running from the tests/ directory
sys.path.append(str(Path(__file__).parent.parent.parent))

from augur import tree

@pytest.fixture
def alignment():
    return Bio.AlignIO.read("tests/data/aa-seq_h3n2_ha_2y_HA1.fasta", "fasta")

@pytest.fixture
def exclude_sites_bed(tmpdir):
    df = pd.DataFrame([
        [42, 45],
        [40, 41],
        [158367030, 158367031],
        [40, 41],
    ])
    filename = str(tmpdir / "exclude_sites.bed")
    df.to_csv(filename, sep='\t')
    return filename

@pytest.fixture
def exclude_sites_txt(tmpdir, mocker):
    filename = str(tmpdir / "exclude_sites.txt")
    with open(filename, "w") as f:
        f.write("618\n")
        f.write("617\n")
        f.write("617\n")
    return filename

@pytest.fixture
def exclude_sites_drm(tmpdir):
    filename = str(tmpdir / "exclude_sites.drm")
    with open(filename, "w") as f:
        f.write("site\tvalue\n")
        f.write("site1\t618\n")
        f.write("site2\t617\n")
        f.write("site3\t618\n")
    return filename

def test_load_excluded_sites_empty_file():
    assert tree.load_excluded_sites(None).tolist() == []

def test_load_excluded_sites_bed(exclude_sites_bed):
    # load_excluded_sites should treat BED files as 0-start,half-open intervals
    # http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
    assert tree.load_excluded_sites(exclude_sites_bed).tolist() == [40, 42, 43, 44, 158367030]

def test_load_excluded_sites_txt(exclude_sites_txt):
    assert tree.load_excluded_sites(exclude_sites_txt).tolist() == [616, 617]

def test_load_excluded_sites_drm(exclude_sites_drm):
    assert tree.load_excluded_sites(exclude_sites_drm).tolist() == [616, 617]
