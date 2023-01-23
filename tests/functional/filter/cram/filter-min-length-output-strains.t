Setup

  $ source "$TESTDIR"/_setup.sh

Filter using only metadata and save results as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-strains filtered_strains.txt > /dev/null

Output should include only the 8 sequences matching the filters (without a header line).

  $ wc -l filtered_strains.txt
  \s*8 .* (re)
