Setup

  $ source "$TESTDIR"/_setup.sh

We don't allow proximity sampling where the query set are themselves a proximal sample

  $ cat >samples1.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >   prox:
  >     query_sample: outbreak
  >   prox_from_prox:
  >     query_sample: prox
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/align.fasta.xz \
  >   --config samples1.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples1.yaml'...
  ERROR: The proximity sample 'prox_from_prox' declares a query sample of 'prox'
  which must itself be a non-proximity sample. I.e. you can't proximity sample directly from another
  proximity sample!
  
  [2]


We don't allow proximity sampling where the context set is a proximal sample

  $ cat >samples1.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Dominican Republic'
  >   prox_one:
  >     query_sample: outbreak
  >   prox_two:
  >     query_sample: outbreak
  >     target_sample: prox_one
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/align.fasta.xz \
  >   --config samples1.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples1.yaml'...
  ERROR: The proximity sample 'prox_two' declares a context sample of 'prox_one'
  which must itself be a non-proximity sample. I.e. you can't proximity sample directly from another
  proximity sample!
  
  [2]

Proximity subsampling requires aligned sequences
(This example is from proximal-subsampling.t but using unaligned sequences)
Because there are multiple outbreak samples, it's the unaligned query seqs which
trigger the error.

  $ cat >samples3.yaml <<~~
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
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples3.yaml \
  >   --output-metadata subsampled3.tsv \
  >   --output-sequences subsampled3.fasta \
  >   --seed 0
  Validating schema of 'samples3.yaml'...
  [outbreak] 10 strains were dropped during filtering
  [outbreak] 	10 were filtered out by the query: "country == 'Dominican Republic'"
  [outbreak] 2 strains passed all filters
  [prox] ERROR: When reading query samples, 'DOM/2016/BB_0183' is 10,621nt, which differs
  [prox] from previous query sequences (10,035nt). Query sequences must be aligned
  [prox] for proximity calculations.
  ERROR: Sample 'prox' failed, see error above.
  [2]

The same as above, but using a single outbreak sequence so the query is "aligned"
(as it's n=1) and thus the error is thrown when reading the --sequences

  $ cat >samples3.yaml <<~~
  > samples:
  >   outbreak:
  >     query: country == 'Puerto Rico'
  >   prox:
  >     query_sample: outbreak
  >     k: 1
  >     max_distance: 1000
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples3.yaml \
  >   --output-metadata subsampled3.tsv \
  >   --output-sequences subsampled3.fasta \
  >   --seed 0
  Validating schema of 'samples3.yaml'...
  [outbreak] 11 strains were dropped during filtering
  [outbreak] 	11 were filtered out by the query: "country == 'Puerto Rico'"
  [outbreak] 1 strain passed all filters
  \[prox\] Read 1 query sequences from .* (re)
  [prox] ERROR: When reading sequences, 'PAN/CDC_259359_V1_V3/2015' has length 10,771nt, which differs
  [prox] from previous query sequences (10,675nt). All sequences must be aligned and
  [prox] the same length for proximity calculations.
  ERROR: Sample 'prox' failed, see error above.
  [2]
