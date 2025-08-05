Setup

  $ source "$TESTDIR"/_setup.sh

Test subsampling with --nthreads 2 for parallel processing.

  $ cat >samples.yaml <<~~
  > samples:
  >   early:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 1
  >     max_date: 2016-02-29
  >     exclude_where:
  >     - region!=South America
  >   recent:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 1
  >     min_date: 2016-03-01
  >     exclude_where:
  >     - region!=South America
  > ~~

Run the subsample command with --nthreads 2.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --subsample-seed 0 \
  >   --nthreads 2
  Validating schema of 'samples.yaml'...
  [early] Running sample 'early'\xe2\x80\xa6 (esc)
  [early] Sampling at 1 per group.
  [early] 11 strains were dropped during filtering
  [early] \t6 were dropped because of 'region!=South America' (esc)
  [early] \t4 were dropped because they were later than 2016.16 or missing a date (esc)
  [early] \t1 was dropped because of subsampling criteria (esc)
  [early] 1 strain passed all filters
  [recent] Running sample 'recent'\xe2\x80\xa6 (esc)
  [recent] Sampling at 1 per group.
  [recent] 11 strains were dropped during filtering
  [recent] \t6 were dropped because of 'region!=South America' (esc)
  [recent] \t3 were dropped because they were earlier than 2016.17 or missing a date (esc)
  [recent] \t2 were dropped because of subsampling criteria (esc)
  [recent] 1 strain passed all filters
  11 strains were dropped during filtering
  \t1 had no metadata (esc)
  \t12 were dropped by `--exclude-all` (esc)
  \\t1 was added back because it was in .*sample_early.* (re)
  \\t1 was added back because it was in .*sample_recent.* (re)
  2 strains passed all filters
