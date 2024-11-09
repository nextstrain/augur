Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree with excluded sites using a compressed input file.

  $ ${AUGUR} tree \
  >  --alignment "$TMP/data/aligned.fasta.xz" \
  >  --exclude-sites "$TMP/data/excluded_sites.txt" \
  >  --output tree_raw.nwk &> /dev/null
