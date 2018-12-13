"""Calculate LBI for a given tree and one or more sets of parameters.
"""
import Bio
import Bio.Phylo
from collections import defaultdict
import json
import numpy as np


def select_nodes_in_season(tree, timepoint, time_window=0.6, **kwargs):
    """Annotate a boolean to each node in the tree if it is alive at the given
    timepoint or prior to the timepoint by the given time window preceding.

    This annotation is used by the LBI and epitope cross-immunity predictors.
    """
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            if node.attr['num_date'] <= timepoint and node.attr['num_date'] > timepoint - time_window:
                node.alive=True
            else:
                node.alive=False
        else:
            node.alive = any(ch.alive for ch in node.clades)


def calculate_LBI(tree, attr="lbi", tau=0.4, transform=lambda x:x, **kwargs):
    '''
    traverses the tree in postorder and preorder to calculate the
    up and downstream tree length exponentially weighted by distance.
    then adds them as LBI
    tree     -- biopython tree for whose node the LBI is being computed
    attr     -- the attribute name used to store the result
    '''
    # Calculate clock length.
    tree.root.clock_length = 0.0
    for node in tree.find_clades():
        for child in node.clades:
            child.clock_length = child.attr['num_date'] - node.attr['num_date']

    # traverse the tree in postorder (children first) to calculate msg to parents
    for node in tree.find_clades(order="postorder"):
        node.down_polarizer = 0
        node.up_polarizer = 0
        for child in node.clades:
            node.up_polarizer += child.up_polarizer
        bl =  node.clock_length / tau
        node.up_polarizer *= np.exp(-bl)
        if node.alive: node.up_polarizer += tau*(1-np.exp(-bl))

    # traverse the tree in preorder (parents first) to calculate msg to children
    for node in tree.get_nonterminals():
        for child1 in node.clades:
            child1.down_polarizer = node.down_polarizer
            for child2 in node.clades:
                if child1!=child2:
                    child1.down_polarizer += child2.up_polarizer

            bl =  child1.clock_length / tau
            child1.down_polarizer *= np.exp(-bl)
            if child1.alive: child1.down_polarizer += tau*(1-np.exp(-bl))

    # go over all nodes and calculate the LBI (can be done in any order)
    max_LBI = 0.0
    for node in tree.find_clades(order="postorder"):
        tmp_LBI = node.down_polarizer
        for child in node.clades:
            tmp_LBI += child.up_polarizer

        node.attr[attr] = transform(tmp_LBI)
        if node.attr[attr] > max_LBI:
            max_LBI = node.attr[attr]

    # Normalize LBI to range [0, 1].
    for node in tree.find_clades():
        node.attr[attr] /= max_LBI
        setattr(node, attr, node.attr[attr])


def register_arguments(parser):
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--branch-lengths", help="JSON with branch lengths and internal node dates estimated by TreeTime", required=True)
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)
    parser.add_argument("--attribute-names", nargs="+", help="names to store distances associated with the corresponding masks", required=True)
    parser.add_argument("--tau", nargs="+", type=float, help="tau value(s) defining the neighborhood of each clade", required=True)
    parser.add_argument("--window", nargs="+", type=float, help="time window(s) to calculate LBI across", required=True)


def run(args):
    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load branch lengths.
    with open(args.branch_lengths, "r") as json_fh:
        branch_lengths = json.load(json_fh)

    # Annotate branch lengths and dates onto tree nodes.
    for node in tree.find_clades():
        node.attr = branch_lengths["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Find maximum time point in the given tree.
    timepoint = max(node.attr["num_date"] for node in tree.find_clades())

    # Calculate LBI for all requested sets of parameters and annotate LBI values
    # to the corresponding attribute names.
    lbi_by_node = defaultdict(dict)
    for i in range(len(args.tau)):
        tau = args.tau[i]
        window = args.window[i]
        attribute_name = args.attribute_names[i]

        # Select nodes that are alive in the given time window.
        select_nodes_in_season(tree, timepoint, window)

        # Calculate LBI.
        calculate_LBI(tree, attribute_name, tau)

        # Collect LBI values into a per-node JSON for export.
        for node in tree.find_clades():
            lbi_by_node[node.name][attribute_name] = node.attr[attribute_name]

    # Export LBI to JSON.
    with open(args.output, "w") as oh:
        json.dump({"nodes": lbi_by_node}, oh, indent=1, sort_keys=True)
