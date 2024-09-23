Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree with an input file that doesn't end in .fasta, and ensure it's not overwritten.

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fa" \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null

  $ sha256sum "$TESTDIR/../data/aligned.fa" | awk '{print $1}'
  169a9f5f70b94e26a2c4ab2b3180d4b463112581438515557a9797adc834863d
