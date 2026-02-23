Setup

  $ source "$TESTDIR"/_setup.sh

Test basic subsampling functionality with a simple config.

  $ cat >samples.yaml <<~~
  > samples:
  >   early:
  >     max_sequences: 1
  >     max_date: 2016-02-29
  >   recent:
  >     max_sequences: 1
  >     min_date: 2016-03-01
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
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

Check that two sequences remain in outputs.

  $ grep -c '^>' subsampled.fasta
  2
  $ tail -n +2 subsampled.tsv | wc -l
  \s*2 (re)
