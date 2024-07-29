Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Build a tree, augmenting existing default arguments with custom arguments.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --tree-builder-args="-czb" \
  >  --output tree_raw.nwk > /dev/null
