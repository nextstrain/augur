Setup

  $ source "$TESTDIR"/_setup.sh

Test config nested under a top level key.

  $ cat >nested_config.yaml <<~~
  > other_tool:
  >   setting: value
  > subsample:
  >   samples:
  >     early:
  >       subsample_max_sequences: 1
  >       max_date: 2016-02-29
  >     recent:
  >       subsample_max_sequences: 1
  >       min_date: 2016-03-01
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-root subsample \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --subsample-seed 0
  Validating schema of 'nested_config.yaml'...
  [early] 11 strains were dropped during filtering
  [early] 	7 were dropped because they were later than 2016.16 or missing a date
  [early] 	4 were dropped because of subsampling criteria
  [early] 1 strain passed all filters
  [recent] 11 strains were dropped during filtering
  [recent] 	5 were dropped because they were earlier than 2016.17 or missing a date
  [recent] 	6 were dropped because of subsampling criteria
  [recent] 1 strain passed all filters
  11 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t1 was added back because it was in .*sample_early.* (re)
  \\t1 was added back because it was in .*sample_recent.* (re)
  2 strains passed all filters

An invalid config root results in an error.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-root invalid \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --subsample-seed 0
  ERROR: Config root key 'invalid' not found in 'nested_config.yaml'
  [2]
