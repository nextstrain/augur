Setup

  $ source "$TESTDIR"/_setup.sh

Defaults are applied to all samples.

  $ cat >config.yaml <<~~
  > defaults:
  >   min_date: 2016
  > samples:
  >   focal:
  >     query: region == 'North America'
  >   contextual:
  >     query: region != 'North America'
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata subsampled.tsv
  Validating schema of 'config.yaml'...
  [focal] 9 strains were dropped during filtering
  [focal] 	8 were filtered out by the query: "region == 'North America'"
  [focal] 	1 was dropped because it was earlier than 2016.0 or missing a date
  [focal] 3 strains passed all filters
  [contextual] 6 strains were dropped during filtering
  [contextual] 	4 were filtered out by the query: "region != 'North America'"
  [contextual] 	2 were dropped because they were earlier than 2016.0 or missing a date
  [contextual] 6 strains passed all filters
  [collect samples] combining outputs from 2 samples
  [collect samples] 3 strains were dropped during filtering
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t3 were added back because they were in .*sample_focal.* (re)
  \[collect samples\] \\t6 were added back because they were in .*sample_contextual.* (re)
  [collect samples] 9 strains passed all filters

Defaults can be overridden by samples.

  $ cat >config.yaml <<~~
  > defaults:
  >   min_date: 2016
  > samples:
  >   focal:
  >     query: region == 'North America'
  >     min_date: 2016-05-01
  >   contextual:
  >     query: region != 'North America'
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata subsampled.tsv
  Validating schema of 'config.yaml'...
  [focal] 11 strains were dropped during filtering
  [focal] 	8 were filtered out by the query: "region == 'North America'"
  [focal] 	3 were dropped because they were earlier than 2016.33 or missing a date
  [focal] 1 strain passed all filters
  [contextual] 6 strains were dropped during filtering
  [contextual] 	4 were filtered out by the query: "region != 'North America'"
  [contextual] 	2 were dropped because they were earlier than 2016.0 or missing a date
  [contextual] 6 strains passed all filters
  [collect samples] combining outputs from 2 samples
  [collect samples] 5 strains were dropped during filtering
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t1 was added back because it was in .*sample_focal.* (re)
  \[collect samples\] \\t6 were added back because they were in .*sample_contextual.* (re)
  [collect samples] 7 strains passed all filters

Some sample options are not supported under defaults.
This is flagged during config validation.

  $ cat >config.yaml <<~~
  > defaults:
  >   min_date: 2016
  >   group_by: year
  > samples:
  >   focal:
  >     query: region == 'North America'
  >     min_date: 2016-05-01
  >   contextual:
  >     query: region != 'North America'
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata subsampled.tsv
  Validating schema of 'config.yaml'...
    .defaults failed: Unexpected property 'group_by'
  ERROR: Validation of 'config.yaml' failed.
  [2]
