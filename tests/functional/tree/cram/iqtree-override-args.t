Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree, replacing existing default arguments with custom arguments.
Since the following custom arguments are incompatible with the default IQ-TREE arguments, this command will only work with the `--override-default-args` flag.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment "$TESTDIR/../data/full_aligned.fasta" \
  >  --tree-builder-args="--polytomy -bb 1000 -bnni" \
  >  --override-default-args \
  >  --output tree_raw.nwk > /dev/null
