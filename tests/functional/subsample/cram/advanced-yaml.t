Setup

  $ source "$TESTDIR"/_setup.sh

YAML anchors and aliases are supported.

Note: this example is an old workaround before adding support for defaults.
Keeping it around as a valid test for advanced YAML syntax.

  $ cat >include.txt <<~~
  > ZKC2/2016
  > VEN/UF_1/2016
  > ~~

  $ cat >config.yaml <<~~
  > subsample_defaults: &subsample_defaults
  >   include:
  >     - include.txt
  > 
  > subsample:
  >   samples:
  >     south_america:
  >       <<: *subsample_defaults
  >       query: region == 'South America'
  >       max_sequences: 1
  >     oceania:
  >       <<: *subsample_defaults
  >       query: region == 'Oceania'
  >       max_sequences: 1
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --config-section subsample \
  >   --output-metadata subsampled.tsv \
  >   --seed 0
  Validating schema of 'config.yaml'...
  [south_america] 9 strains were dropped during filtering
  [south_america] 	6 were filtered out by the query: "region == 'South America'"
  \[south_america\] \\t2 were added back because they were in .*/include.txt.* (re)
  [south_america] 	5 were dropped because of subsampling criteria
  [south_america] 3 strains passed all filters
  [oceania] 10 strains were dropped during filtering
  [oceania] 	11 were filtered out by the query: "region == 'Oceania'"
  \[oceania\] \\t2 were added back because they were in .*/include.txt.* (re)
  [oceania] 	0 were dropped because of subsampling criteria
  [oceania] 2 strains passed all filters
  [collect samples] combining outputs from 2 samples
  [collect samples] 9 strains were dropped during filtering
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t3 were added back because they were in .*south_america.* (re)
  \[collect samples\] \\t2 were added back because they were in .*oceania.* (re)
  [collect samples] 3 strains passed all filters

The ids in include.txt are present in the output.

  $ grep -c 'ZKC2/2016' subsampled.tsv
  1
  $ grep -c 'VEN/UF_1/2016' subsampled.tsv
  1
