Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter using only metadata without a sequence index.
This should work because the requested filters don't rely on sequence information.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ rm -f "$TMP/filtered_strains.txt"
