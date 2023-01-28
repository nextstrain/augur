Setup

  $ source "$TESTDIR"/_setup.sh

Filter using only metadata without a sequence index.
This should work because the requested filters don't rely on sequence information.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --output-strains filtered_strains.txt > /dev/null
