Setup

  $ source "$TESTDIR"/_setup.sh

Try building a tree with IQ-TREE with more threads (4) than there are input sequences (3).

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 4 > /dev/null
