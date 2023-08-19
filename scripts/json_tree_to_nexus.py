from __future__ import print_function
import sys
sys.path.append('..')
import json
from base.io_util import json_to_tree
from base.utils import save_as_nexus
import Bio.Phylo
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    parser.add_argument('-t', '--tree', required=True, type=str, help="JSON tree file")
    parser.add_argument('--resolve', action="store_true", help="resolve polytomies")
    parser.add_argument('--temporal', action="store_true", help="export branch lengths in time")
    parser.add_argument('-o', '--outfile', type=str, default=None, help="output nexus filename (optional)")

    return parser.parse_args()


def resolve_polytomies(tree):

    def assign_attrs(node, source):
        # Assign all non-children attributes.
        for attr, value in source.__dict__.items():
            if attr != "children":
                setattr(node, attr, value)
        # node.xvalue = 0
        # node.attr["div"] = 0

    def order_children_in_place(clades):
        clades.sort(key=lambda x: x.attr["div"])
        clades.sort(key=lambda x: x.attr["num_date"])

    def ladderise_node(node):
        """input: node object with 2 or more chilred
        returns node object but with 2 children. One of which may be a polytomy"""
        if len(node.clades) == 2:
            return node
        order_children_in_place(node.clades)
        new_node = Bio.Phylo.Newick.Clade()
        if "resolved" in node.name:
            tmp = node.name.split("_")
            new_node.name = "_".join(tmp[:-1]) + "_" + str(int(tmp[-1] + 1))
        else:
            new_node.name = node.strain + "_resolved_1"
        assign_attrs(new_node, node) # modifies in place
        new_node.clades = node.clades[1:]
        node.clades = [node.clades[0], ladderise_node(new_node)]
        return node

    def traverse_and_resolve_first(tree):
        for clade in tree.find_clades(terminal=False, order="preorder"):
            if len(clade.clades) != 2:
                print("resolving ", clade.name, " with ", len(clade.clades), " children using temporal data")
                clade = ladderise_node(clade)
                return True
        return False

    ## fn start
    polytomies_exist = True; # stopping condition
    while polytomies_exist:
        polytomies_exist = traverse_and_resolve_first(tree)
    assert(tree.is_bifurcating())
    tree.ladderize()

if __name__=="__main__":
    params = collect_args()
    outfile = params.outfile if params.outfile != None else params.tree.split(".json")[0]+".nex"

    json_fh = open(params.tree, "r")
    json_dict = json.load(json_fh)
    tree = json_to_tree(json_dict)
    if params.resolve:
        resolve_polytomies(tree)
    if params.temporal:
        save_as_nexus(tree, outfile, "num_date")
    else:
        save_as_nexus(tree, outfile, "div")
