Setup

  $ source "$TESTDIR"/_setup.sh

Try to group data without any grouping arguments.
This should fail with a helpful error message.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by year month \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  [2]
