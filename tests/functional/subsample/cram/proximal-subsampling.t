Setup

  $ source "$TESTDIR"/_setup.sh

Test basic proximal subsampling: an outbreak sample filtered to Dominican Republic,
then a proximity sample finding k=1 closest sequences using the outbreak as focal.

  $ cat >samples1.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >   prox:
  >     focal_sample: outbreak
  >     k: 1
  >     max_distance: 1000
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/align.fasta.xz \
  >   --config samples1.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples1.yaml'...
  [outbreak] 10 strains were dropped during filtering
  [outbreak] \t10 were filtered out by the query: "country == 'Dominican Republic'" (esc)
  [outbreak] 2 strains passed all filters
  \[prox\] Read 2 focal sequences from .* (re)
  \[prox\] Loaded 10 comparison sequences from .*. Excluded 2 as they were in the focal sequences. (re)
  \[prox\]\s* (re)
  \[prox\] Proximity calculations complete. Total time: .*s (re)
  [prox] Found 2 unique neighbor strains for 2 focal strains
  \[prox\] Wrote 2 strains to .* (re)
  [collect samples] combining outputs from 2 samples
  [collect samples] 9 strains were dropped during filtering
  [collect samples] \t1 had no metadata (esc)
  [collect samples] \t12 were dropped by `--exclude-all` (esc)
  \[collect samples\] \\t2 were added back because they were in .*sample_outbreak.* (re)
  \[collect samples\] \\t2 were added back because they were in .*sample_prox.* (re)
  [collect samples] 4 strains passed all filters

Test proximal subsampling with explicit context sample: an outbreak sample filtered
to Dominican Republic, a context sample filtered to South America, and a proximity
sample finding k=1 closest from the context using the outbreak as focal.
The expectation is a final output with 2 from North America and 2 from South America

  $ cat >samples2.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >   ctx:
  >     query: region == 'South America'
  >     drop_sample: true
  >   prox:
  >     focal_sample: outbreak
  >     context_sample: ctx
  >     k: 1
  >     max_distance: 1000
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/align.fasta.xz \
  >   --config samples2.yaml \
  >   --output-metadata subsampled2.tsv \
  >   --output-sequences subsampled2.fasta \
  >   --seed 0
  Validating schema of 'samples2.yaml'...
  [outbreak] 10 strains were dropped during filtering
  [outbreak] \t10 were filtered out by the query: "country == 'Dominican Republic'" (esc)
  [outbreak] 2 strains passed all filters
  [ctx] 6 strains were dropped during filtering
  [ctx] \t6 were filtered out by the query: "region == 'South America'" (esc)
  [ctx] 6 strains passed all filters
  \[prox\] Read 2 focal sequences from .* (re)
  \[prox\] Loaded 6 comparison sequences from .*\. Excluded 0 as they were in the focal sequences\. (re)
  \[prox\]\s* (re)
  \[prox\] Proximity calculations complete. Total time: .*s (re)
  [prox] Found 2 unique neighbor strains for 2 focal strains
  \[prox\] Wrote 2 strains to .* (re)
  [collect samples] combining outputs from 2 samples
  [collect samples] 9 strains were dropped during filtering
  [collect samples] \t1 had no metadata (esc)
  [collect samples] \t12 were dropped by `--exclude-all` (esc)
  \[collect samples\] \\t2 were added back because they were in .*sample_outbreak.* (re)
  \[collect samples\] \\t2 were added back because they were in .*sample_prox.* (re)
  [collect samples] 4 strains passed all filters

  $ tail -n +2 subsampled2.tsv | cut -f 5 | sort | uniq -c | sed 's/^ */  /'
    2 North America
    2 South America


We can subsample further from a proximal set. The use cases for this are (I think)
rare, but with future improvements such as exporting proximity priorities it could
become more common.

  $ cat >samples3.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >     drop_sample: true
  >   prox:
  >     focal_sample: outbreak
  >     k: 5
  >     max_distance: 1000
  >     drop_sample: true
  >   prox-in-north-america:
  >     context_sample: prox
  >     query: region == 'North America'
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/align.fasta.xz \
  >   --config samples3.yaml \
  >   --output-metadata subsampled3.tsv \
  >   --output-sequences subsampled3.fasta \
  >   --seed 0 &> /dev/null

  $ tail -n +2 subsampled3.tsv | cut -f 1
  PRVABC59
