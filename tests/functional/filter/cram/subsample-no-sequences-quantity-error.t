Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to group data without any grouping arguments.
This should fail with a helpful error message.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by year month \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  [2]

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by year month \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  [2]
