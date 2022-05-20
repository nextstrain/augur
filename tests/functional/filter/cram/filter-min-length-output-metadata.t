Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter using only metadata without sequence input or output and save results as filtered metadata.

  $ ${AUGUR} filter \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null

Output should include the 8 sequences matching the filters and a header line.

  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*9 .* (re)
  $ rm -f "$TMP/filtered_metadata.tsv"
