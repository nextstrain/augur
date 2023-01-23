Setup

  $ source "$TESTDIR"/_setup.sh

Filter using only metadata without sequence input or output and save results as filtered metadata.

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-metadata filtered_metadata.tsv > /dev/null

Output should include the 8 sequences matching the filters and a header line.

  $ wc -l filtered_metadata.tsv
  \s*9 .* (re)
