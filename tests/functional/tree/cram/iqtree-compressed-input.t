Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree with excluded sites using a compressed input file.

  $ cp "$TESTDIR/../data/aligned.fasta.xz" .
  $ cp "$TESTDIR/../data/excluded_sites.txt" .
  $ ${AUGUR} tree \
  >  --alignment aligned.fasta.xz \
  >  --exclude-sites excluded_sites.txt \
  >  --output tree_raw.nwk &> /dev/null
