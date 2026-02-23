Setup

  $ source "$TESTDIR"/_setup.sh

Test config nested under a top level key.

  $ cat >nested_config.yaml <<~~
  > other_tool:
  >   setting: value
  > subsample:
  >   samples:
  >     early:
  >       max_sequences: 1
  >       max_date: 2016-02-29
  >     recent:
  >       max_sequences: 1
  >       min_date: 2016-03-01
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-section subsample \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --seed 0
  Validating schema of 'nested_config.yaml'...
  [early] 11 strains were dropped during filtering
  [early] 	7 were dropped because they were later than 2016.16 or missing a date
  [early] 	4 were dropped because of subsampling criteria
  [early] 1 strain passed all filters
  [recent] 11 strains were dropped during filtering
  [recent] 	5 were dropped because they were earlier than 2016.17 or missing a date
  [recent] 	6 were dropped because of subsampling criteria
  [recent] 1 strain passed all filters
  [collect samples] combining outputs from 2 samples
  [collect samples] 11 strains were dropped during filtering
  [collect samples] 	1 had no metadata
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t1 was added back because it was in .*sample_early.* (re)
  \[collect samples\] \\t1 was added back because it was in .*sample_recent.* (re)
  [collect samples] 2 strains passed all filters

An invalid config path results in an error.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-section invalid \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --seed 0
  ERROR: Config section 'invalid' not found in 'nested_config.yaml'
  [2]

Test config nested under multiple keys.

  $ cat >nested_config.yaml <<~~
  > builds:
  >   avian-flu/h5n1/ha/all-time:
  >     other_tool:
  >       setting: value
  >     subsample:
  >       samples:
  >         early:
  >           max_sequences: 1
  >           max_date: 2016-01-01
  >         recent:
  >           max_sequences: 1
  >           min_date: 2016
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-section builds avian-flu/h5n1/ha/all-time subsample \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --seed 0
  Validating schema of 'nested_config.yaml'...
  [early] 11 strains were dropped during filtering
  [early] 	9 were dropped because they were later than 2016.0 or missing a date
  [early] 	2 were dropped because of subsampling criteria
  [early] 1 strain passed all filters
  [recent] 11 strains were dropped during filtering
  [recent] 	3 were dropped because they were earlier than 2016.0 or missing a date
  [recent] 	8 were dropped because of subsampling criteria
  [recent] 1 strain passed all filters
  [collect samples] combining outputs from 2 samples
  [collect samples] 11 strains were dropped during filtering
  [collect samples] 	1 had no metadata
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t1 was added back because it was in .*sample_early.* (re)
  \[collect samples\] \\t1 was added back because it was in .*sample_recent.* (re)
  [collect samples] 2 strains passed all filters

An invalid config path results in an error.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-section builds avian-flu/h5n1/ha/6m subsample \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --seed 0
  ERROR: Config section 'builds' → 'avian-flu/h5n1/ha/6m' not found in 'nested_config.yaml'
  [2]
