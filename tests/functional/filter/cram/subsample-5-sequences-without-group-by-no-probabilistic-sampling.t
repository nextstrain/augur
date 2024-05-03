Setup

  $ source "$TESTDIR"/_setup.sh

Filter with subsampling where no more than 5 sequences are requested and no groups are specified.
This generates a dummy category and subsamples from there. With no-probabilistic-sampling we expect exactly 5 sequences.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output filtered.fasta 2>/dev/null
  $ grep ">" filtered.fasta | wc -l
  \s*5 (re)
