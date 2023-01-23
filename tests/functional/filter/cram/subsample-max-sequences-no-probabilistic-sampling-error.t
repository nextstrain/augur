Setup

  $ source "$TESTDIR"/_setup.sh

This should fail, as probabilistic sampling is explicitly disabled.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output filtered.fasta
  ERROR: Asked to provide at most 5 sequences, but there are 8 groups.
  [2]
