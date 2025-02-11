import sys
from augur.utils import read_tree

if __name__ == "__main__":
    T = read_tree(sys.argv[1])
    root_child_tips = [c for c in T.root.clades if c.is_terminal()]
    if len(root_child_tips)==0:
        print("No children of the root were terminal nodes")
    elif len(root_child_tips)>1:
        print("Multiple children of the root were terminal nodes")
    else:
        print(f"Tree root has a single terminal child {root_child_tips[0].name!r}")
