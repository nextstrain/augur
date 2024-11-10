Setup

  $ source "$TESTDIR"/_setup.sh

Try building a tree with IQ-TREE using its ModelTest functionality, by supplying a substitution model of "auto".

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --method iqtree \
  >  --substitution-model auto \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null
