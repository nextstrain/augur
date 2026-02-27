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
  >     context_sample: prox_one
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
