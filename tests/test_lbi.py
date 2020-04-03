"""
Unit tests for LBI calculation
"""

import Bio
import Bio.Phylo
from collections import defaultdict
import json
import numpy as np
import pytest

import os
from os import path

from augur.lbi import select_nodes_in_season, calculate_LBI, register_arguments, run
from augur.utils import json_to_tree, write_json


@pytest.fixture
def tree():
    """Returns an annotated Bio.Phylo tree.
    """
    with open("data/flu_seasonal_h3n2_ha_3y_tree.json") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)
    return tree


# Define global variables for test_select_nodes_in_season
TIMEPOINT = 10
TIME_WINDOW = 0.6

def test_select_nodes_in_season(tree, timepoint=TIMEPOINT, time_window=TIME_WINDOW):
    """Test boolean annotations of tree.
    """

    select_nodes_in_season(tree, timepoint=TIMEPOINT, time_window=TIME_WINDOW)

    clades = list(tree.find_clades(order="postorder"))

    
    # Make sure that tree contains nodes
    assert len(clades) > 0

    # Check to make sure all nodes have 'num_date' attribute
    assert all(node.attr['num_date'] for node in tree.find_clades(order="postorder"))

    # Check to make sure all nodes have an alive attr
    assert all(hasattr(node, "alive") for node in tree.find_clades(order="postorder"))


## Test the calculation of LBI

# Define global variables for test_calculate_LBI 

ATTR = "lbi"
TAU = 0.4
TRANSFORM = lambda x:x
NORMALIZE = True

def test_calculate_LBI(tree, attr=ATTR, tau=TAU, transform=TRANSFORM, normalize=NORMALIZE):
    """Test calculation of LBI.
    """

    # Run this function to add alive attr to nodes
    select_nodes_in_season(tree, timepoint=TIMEPOINT, time_window=TIME_WINDOW)

    # Run the calculate_LBI function
    calculate_LBI(tree, attr=ATTR, tau=TAU, transform=TRANSFORM, normalize=NORMALIZE)

    # Check to make sure all nodes alive attr is either True of False
    assert all(hasattr(node, "attr") for node in tree.find_clades())

    # Check that the LBI has been normalized between 0 and 1
    assert all([node.attr for node in tree.find_clades()]) >= 0 and all([node.attr for node in tree.find_clades()]) <= 1


### Test the run function
## Define necessary variables for test_run function
ATTRIBUTE_NAMES = list(str(e) for e in np.arange(0,35))

# Set up lists of tau and window values equal to length of Zika fasta
TAU_LIST = []
for i in range(0,35):
    TAU_LIST.append(0.4)

WINDOWS = []
for i in range(0,35):
    WINDOWS.append(0.6)

# Construct a mock Args class to mimic the parser input for run function
class Args():
    """Mock args class for input to run function."""

    def __init__(self):

        self.tree = "data/zika-tutorial/results/tree.nwk"
        self.branch_lengths = "data/zika-tutorial/results/branch_lengths.json"
        self.tau = TAU_LIST
        self.window = WINDOWS
        self.attribute_names = ATTRIBUTE_NAMES
        self.no_normalization = True
        self.output = "data/zika-tutorial/results/lbi_test_run_output"


@pytest.fixture
def args():
    """Instantiate an Args object"""
    args = Args()

    return args


def test_run(args):
    """Test the run function.
    """

    # Run the run function with test args
    run(args)

    # Check for the final output file
    assert path.exists("data/zika-tutorial/results/lbi_test_run_output")

    # Check to make sure the output file is in the correct format
    json_file = open("data/zika-tutorial/results/lbi_test_run_output")
    json_dict = json.load(json_file)

    # Check that the json file has the correct header
    assert "generated_by" in json_dict.keys()

    # Remove the test_run output file
    os.remove("data/zika-tutorial/results/lbi_test_run_output")
    
















