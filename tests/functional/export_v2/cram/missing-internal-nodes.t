Setup

  $ source "$TESTDIR"/_setup.sh

Test a tree with unlabelled internal nodes.
In Augur v24.3.0 and earlier these would work _sometimes_ depending
on what code paths were hit, but in most real-life cases would raise
an uncaught error. We now check for these when we parse the tree.

  $ echo "(tipA:1,(tipB:1,tipC:1):2,(tipD:3,tipE:4,tipF:1):5):0;" > tree1.nwk

  $ ${AUGUR} export v2 \
  >   --tree tree1.nwk \
  >   --metadata "$TESTDIR/../data/dataset1_metadata_with_name.tsv" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" \
  >   --output test1.json
  ERROR: Tree contains unnamed nodes.+ (re)
  [2]
