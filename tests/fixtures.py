"""
Fixtures shared by multiple test modules.
"""
import Bio.Align.AlignInfo
import Bio.Phylo
import Bio.SeqIO
import json
import numpy as np
import os
import pytest
import sys

# we assume (and assert) that this script is running from the tests/ directory
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from base.io_util import json_to_tree


@pytest.fixture
def multiple_sequence_alignment():
    """Returns a multiple sequence alignment containing a small test set of H3N2
    sequences.
    """
    msa = Bio.AlignIO.read("tests/data/fitness_model/H3N2_alignment.cleaned.fasta", "fasta")
    return msa


@pytest.fixture
def real_tree(multiple_sequence_alignment):
    """Returns a tree built with FastTree from a small set of nucleotide sequences
    for H3N2.
    """
    # Load the tree.
    tree = Bio.Phylo.read("tests/data/fitness_model/H3N2_tree.newick", "newick")

    # Make a lookup table of name to sequence.
    sequences_by_name = dict([(alignment.name, str(alignment.seq))
                  for alignment in multiple_sequence_alignment])

    # Assign sequences to the tree.
    index = 0
    for node in tree.find_clades():
        if node.name is not None:
            node.sequence = np.fromstring(sequences_by_name[node.name], "S1")

            # Since sequence names look like "A/Singapore/TT0495/2017",
            # convert the last element to a floating point value for
            # simplicity.
            node.attr = {"num_date": float(node.name.split("/")[-1])}
        else:
            # Build a "dumb" consensus from the alignment for the
            # ancestral node and assign an arbitrary date in the
            # past.
            summary = Bio.Align.AlignInfo.SummaryInfo(multiple_sequence_alignment)
            node.sequence = np.fromstring(str(summary.dumb_consensus(threshold=0.5, ambiguous="N")), "S1")
            node.attr = {"num_date": 2014.8}

        node.clade = index
        index += 1

    return tree


@pytest.fixture
def json_tree():
    """Returns an annotated Bio.Phylo tree.
    """
    with open("tests/data/flu_seasonal_h3n2_ha_3y_tree.json") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)
    return tree
