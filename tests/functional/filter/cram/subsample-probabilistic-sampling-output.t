Setup

  $ source "$TESTDIR"/_setup.sh

Check output of probabilistic sampling.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by region year month \
  >  --subsample-max-sequences 3 \
  >  --probabilistic-sampling \
  >  --subsample-seed 314159 \
  >  --output-metadata filtered_metadata.tsv
  WARNING: Asked to provide at most 3 sequences, but there are 8 groups.
  Sampling probabilistically at 0.3750 sequences per group, meaning it is possible to have more than the requested maximum of 3 sequences after filtering.
  10 strains were dropped during filtering
  	1 was dropped during grouping due to ambiguous year information
  	1 was dropped during grouping due to ambiguous month information
  	8 were dropped because of subsampling criteria, using seed 314159
  2 strains passed all filters
