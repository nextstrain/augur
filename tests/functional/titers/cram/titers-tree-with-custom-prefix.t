Setup

  $ source "$TESTDIR"/_setup.sh

Test titer tree model with a custom prefix for the node data attributes in the output.

  $ ${AUGUR} titers tree \
  >   --tree ../data/tree.nwk \
  >   --titers ../data/titers.tsv \
  >   --attribute-prefix custom_prefix_ \
  >   --output $TMP/titers-tree.json > /dev/null
  Read titers from ../data/titers.tsv, found:
   --- 62 strains
   --- 15 data sources
   --- 272 total measurements
  $ grep custom_prefix_cTiter $TMP/titers-tree.json | wc -l
  \s*120 (re)
