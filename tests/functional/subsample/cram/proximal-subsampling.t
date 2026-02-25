Setup

  $ source "$TESTDIR"/_setup.sh

Test basic proximal subsampling: an outbreak sample filtered to Dominican Republic,
then a proximity sample finding k=1 closest sequences using the outbreak as query.

  $ cat >samples1.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >   prox:
  >     query_sample: outbreak
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
  \[prox\] Read 2 query sequences from .* (re)
  \[prox\] Loaded 10 comparison sequences from .*. Excluded 2 as they were in the focal sequences. (re)
  \[prox\]\s* (re)
  \[prox\] Proximity calculations complete. Total time: .*s (re)
  [prox] Found 2 unique neighbor strains for 2 query strains
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
sample finding k=1 closest from the context using the outbreak as query.
The expectation is a final output with 2 from North America and 2 from South America

  $ cat >samples2.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >   ctx:
  >     query: region == 'South America'
  >     drop_sample: true
  >   prox:
  >     query_sample: outbreak
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
  \[prox\] Read 2 query sequences from .* (re)
  \[prox\] Loaded 6 comparison sequences from .*\. Excluded 0 as they were in the focal sequences\. (re)
  \[prox\]\s* (re)
  \[prox\] Proximity calculations complete. Total time: .*s (re)
  [prox] Found 2 unique neighbor strains for 2 query strains
  \[prox\] Wrote 2 strains to .* (re)
  [collect samples] combining outputs from 2 samples
  [collect samples] 9 strains were dropped during filtering
  [collect samples] \t1 had no metadata (esc)
  [collect samples] \t12 were dropped by `--exclude-all` (esc)
  \[collect samples\] \\t2 were added back because they were in .*sample_outbreak.* (re)
  \[collect samples\] \\t2 were added back because they were in .*sample_prox.* (re)
  [collect samples] 4 strains passed all filters

  $ cat subsampled2.tsv | csvtk cut -t --delete-header -f region | sort | uniq -c
     2 North America
     2 South America
