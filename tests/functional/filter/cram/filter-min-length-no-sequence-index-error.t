Setup

  $ source "$TESTDIR"/_setup.sh

Try to filter using only metadata without a sequence index.
This should fail because the requested filters rely on sequence information.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 10000 \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: You need to provide a sequence index or sequences to filter on sequence-specific information.
  [2]
