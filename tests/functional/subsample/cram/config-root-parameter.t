Setup

  $ source "$TESTDIR"/_setup.sh

Test config with --config-root parameter.

  $ cat >nested_config.yaml <<~~
  > other_tool:
  >   setting: value
  > subsample:
  >   samples:
  >     early:
  >       group_by:
  >       - region
  >       subsample_max_sequences: 1
  >       max_date: 2016-02-29
  >       exclude_where:
  >       - region!=South America
  >     recent:
  >       group_by:
  >       - region
  >       subsample_max_sequences: 1
  >       min_date: 2016-03-01
  >       exclude_where:
  >       - region!=South America
  > ~~

Run subsample with --config-root parameter.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_config.yaml \
  >   --config-root subsample \
  >   --output-metadata nested_subsampled.tsv \
  >   --output-sequences nested_subsampled.fasta \
  >   --subsample-seed 0
  Validating schema of 'nested_config.yaml'...
  Sampling at 1 per group.
  11 strains were dropped during filtering
  	6 were dropped because of 'region!=South America'
  	4 were dropped because they were later than 2016.16 or missing a date
  	1 was dropped because of subsampling criteria
  1 strain passed all filters
  Sampling at 1 per group.
  11 strains were dropped during filtering
  	6 were dropped because of 'region!=South America'
  	3 were dropped because they were earlier than 2016.17 or missing a date
  	2 were dropped because of subsampling criteria
  1 strain passed all filters
  11 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t1 was added back because it was in .*sample_early.* (re)
  \\t1 was added back because it was in .*sample_recent.* (re)
  2 strains passed all filters

Check that --config-root produces the same results.

  $ grep -c '^>' nested_subsampled.fasta
  2
  $ tail -n +2 nested_subsampled.tsv | wc -l | tr -d ' '
  2
  $ cut -f5 nested_subsampled.tsv | tail -n +2 | sort | uniq
  South America
