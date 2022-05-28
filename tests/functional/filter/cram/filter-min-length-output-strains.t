Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter using only metadata and save results as a list of filtered strains.
Output should include only the 8 sequences matching the filters (without a header line).

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*8 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*8 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"
