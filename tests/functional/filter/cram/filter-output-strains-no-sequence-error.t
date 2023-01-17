Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to filter with sequence outputs and no sequence inputs.
This should fail.

  $ ${AUGUR} filter \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-length 10000 \
  >  --output "$TMP/filtered.fasta" > /dev/null
  ERROR: You need to provide sequences to output sequences.
  [2]
