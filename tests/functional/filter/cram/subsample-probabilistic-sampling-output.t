Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Check output of probabilistic sampling.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by region year month \
  >  --subsample-max-sequences 3 \
  >  --probabilistic-sampling \
  >  --subsample-seed 314159 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  WARNING: Asked to provide at most 3 sequences, but there are 8 groups.
  Sampling probabilistically at 0.3633 sequences per group, meaning it is possible to have more than the requested maximum of 3 sequences after filtering.
  10 strains were dropped during filtering
  \t1 were dropped during grouping due to ambiguous year information (esc)
  \t1 were dropped during grouping due to ambiguous month information (esc)
  \t8 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  2 strains passed all filters

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by region year month \
  >  --subsample-max-sequences 3 \
  >  --probabilistic-sampling \
  >  --subsample-seed 314159 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  WARNING: Asked to provide at most 3 sequences, but there are 8 groups.
  Sampling probabilistically at 0.3633 sequences per group, meaning it is possible to have more than the requested maximum of 3 sequences after filtering.
  10 strains were dropped during filtering
  \t1 were dropped during grouping due to ambiguous month information (esc)
  \t1 were dropped during grouping due to ambiguous year information (esc)
  \t8 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  2 strains passed all filters
