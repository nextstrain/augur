Setup

  $ source "$TESTDIR"/_setup.sh

Explicitly use probabilistic subsampling to handle the case when there are more available groups than requested sequences.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --probabilistic-sampling \
  >  --output-strains filtered_strains_probabilistic.txt > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 8 groups.
  Sampling probabilistically at 0.6250 sequences per group, meaning it is possible to have more than the requested maximum of 5 sequences after filtering.
  10 strains were dropped during filtering
  	1 had no metadata
  	1 had no sequence data
  	1 was dropped because it was earlier than 2012.0 or missing a date
  	1 was dropped during grouping due to ambiguous month information
  	6 were dropped because of subsampling criteria, using seed 314159
  3 strains passed all filters

Using the default probabilistic subsampling, should work the same as the previous case.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --output-strains filtered_strains_default.txt > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 8 groups.
  Sampling probabilistically at 0.6250 sequences per group, meaning it is possible to have more than the requested maximum of 5 sequences after filtering.
  10 strains were dropped during filtering
  	1 had no metadata
  	1 had no sequence data
  	1 was dropped because it was earlier than 2012.0 or missing a date
  	1 was dropped during grouping due to ambiguous month information
  	6 were dropped because of subsampling criteria, using seed 314159
  3 strains passed all filters

By setting the subsample seed above, we should get the same results for both runs.

  $ diff -u <(sort filtered_strains_probabilistic.txt) <(sort filtered_strains_default.txt)
