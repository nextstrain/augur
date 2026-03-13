Setup

  $ source "$TESTDIR"/_setup.sh

Test subsampling functionality where samples depend on other samples
(This config could be written in a single sample, but we're deliberately using three here)

  $ cat >samples.yaml <<~~
  > samples:
  >   A:
  >     query: region == 'North America'
  >     drop_sample: true
  >   B:
  >     target_sample: A
  >     query: country == 'Dominican Republic'
  >     drop_sample: true
  >   C:
  >     target_sample: B
  >     min_date: 2016-04-10
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  [A] 8 strains were dropped during filtering
  [A] 	8 were filtered out by the query: "region == 'North America'"
  [A] 4 strains passed all filters
  [B] 2 strains were dropped during filtering
  [B] 	2 were filtered out by the query: "country == 'Dominican Republic'"
  [B] 2 strains passed all filters
  [C] 1 strain was dropped during filtering
  [C] 	1 was dropped because it was earlier than 2016.27 or missing a date
  [C] 1 strain passed all filters
  [collect samples] combining outputs from 1 samples
  [collect samples] 12 strains were dropped during filtering
  [collect samples] 	1 had no metadata
  [collect samples] 	12 were dropped by `--exclude-all`
  \[collect samples\] \\t1 was added back because it was in .*sample_C.* (re)
  [collect samples] 1 strain passed all filters
