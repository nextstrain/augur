Setup

  $ source "$TESTDIR"/_setup.sh

Try building a tree with IQ-TREE.

  $ cp "$TESTDIR/../data/aligned.fasta" .
  $ ${AUGUR} tree \
  >  --alignment aligned.fasta \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null
