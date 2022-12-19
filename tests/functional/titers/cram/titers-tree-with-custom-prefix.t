Setup

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../../../bin/augur}"

Test titer tree model with a custom prefix for the node data attributes in the output.

  $ ${AUGUR} titers tree \
  >   --tree ../data/tree.nwk \
  >   --titers ../data/titers.tsv \
  >   --attribute-prefix custom_prefix_ \
  >   --output $TMP/titers-tree.json > /dev/null
  Read titers from ../data/titers.tsv, found:
   --- 61 strains
   --- 15 data sources
   --- 232 total measurements
  $ grep custom_prefix_cTiter $TMP/titers-tree.json | wc -l
  \s*120 (re)
