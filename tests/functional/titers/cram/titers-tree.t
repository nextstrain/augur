Setup

  $ source "$TESTDIR"/_setup.sh

Test titer tree model.

  $ ${AUGUR} titers tree \
  >   --tree ../data/tree.nwk \
  >   --titers ../data/titers.tsv \
  >   --output $TMP/titers-tree.json > /dev/null
  Read titers from ../data/titers.tsv, found:
   --- 62 strains
   --- 15 data sources
   --- 272 total measurements
  $ grep cTiter $TMP/titers-tree.json | wc -l
  \s*120 (re)
