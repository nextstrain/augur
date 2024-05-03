Setup

  $ source "$TESTDIR"/_setup.sh

Filter with subsampling, requesting no more than 8 sequences.
With 8 groups to subsample from (after filtering), this should produce one sequence per group.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 8 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output filtered.fasta 2>/dev/null
  $ grep ">" filtered.fasta | wc -l
  \s*8 (re)
