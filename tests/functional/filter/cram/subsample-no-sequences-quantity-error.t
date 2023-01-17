Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to group data without any grouping arguments.
This should fail with a helpful error message.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by year month \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  [2]
