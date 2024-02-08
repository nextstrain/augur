Setup

  $ source "$TESTDIR"/_setup.sh

Ensure probabilistic sampling is not used when unnecessary.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by region year month \
  >  --subsample-max-sequences 10 \
  >  --probabilistic-sampling \
  >  --subsample-seed 314159 \
  >  --output-metadata filtered_metadata.tsv
  Sampling at 10 per group.
  2 strains were dropped during filtering
  	1 was dropped during grouping due to ambiguous year information
  	1 was dropped during grouping due to ambiguous month information
  	0 were dropped because of subsampling criteria, using seed 314159
  10 strains passed all filters
