  $ source "$TESTDIR"/_setup.sh

  $ echo "(tipA:1,(tipA:1,tipA:1)internalBC:2,(tipD:3,tipE:4,tipF:1)internalDEF:5)ROOT:0;" > tree1.nwk

  $ ${AUGUR} export v2 \
  >   --tree tree1.nwk \
  >   --output minimal.json
  ERROR: 1 node names occur multiple times in the tree: 'tipA'
  [2]

