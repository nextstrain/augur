Setup

  $ source "$TESTDIR"/_setup.sh

Try to filter with sequence outputs and no sequence inputs.
This should fail.

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 10000 \
  >  --output filtered.fasta > /dev/null
  ERROR: You need to provide sequences to output sequences.
  [2]
