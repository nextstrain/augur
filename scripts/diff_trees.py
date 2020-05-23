import argparse
import Bio.Phylo
import deepdiff


def clade_to_items(clade, attrs=("name", "branch_length")):
    """Recursively convert a clade of a tree to a list of nested lists according to
    the topology of the clade with the requested attributes per node.

    >>> from io import StringIO
    >>> treedata = "(A, (B, C), (D, E))"
    >>> handle = StringIO(treedata)
    >>> tree = Bio.Phylo.read(handle, "newick")
    >>> clade_to_items(tree.root)
    [[None, None], [['A', None]], [[None, None], [['B', None]], [['C', None]]], [[None, None], [['D', None]], [['E', None]]]]
    """
    items = [[
        getattr(clade, attr)
        for attr in attrs
    ]]

    for child in clade.clades:
        items.extend([clade_to_items(child)])

    return items


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("first_tree", help="first Newick tree to compare")
    parser.add_argument("second_tree", help="second Newick tree to compare")
    parser.add_argument("--attributes", nargs="+", default=["name", "branch_length"], help="node attributes to include in comparison")
    parser.add_argument("--significant-digits", type=int, default=5, help="number of significant digits to use when comparing branch lengths")

    args = parser.parse_args()

    first_tree = Bio.Phylo.read(args.first_tree, "newick")
    second_tree = Bio.Phylo.read(args.second_tree, "newick")

    first_tree_items = clade_to_items(first_tree.root, attrs=args.attributes)
    second_tree_items = clade_to_items(second_tree.root, attrs=args.attributes)

    print(
        deepdiff.DeepDiff(
            first_tree_items,
            second_tree_items,
            significant_digits=args.significant_digits
        )
    )
