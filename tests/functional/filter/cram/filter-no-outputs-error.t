Setup

  $ source "$TESTDIR"/_setup.sh

Try to filter without any outputs.

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 10000 > /dev/null
  ERROR: You need to select at least one output.
  [2]
