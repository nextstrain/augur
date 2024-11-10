Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree, augmenting existing default arguments with custom arguments.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --tree-builder-args="--polytomy" \
  >  --output tree_raw.nwk > /dev/null
