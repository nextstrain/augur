Setup

  $ source "$TESTDIR"/_setup.sh

A target sample which doesn't exist

  $ cat >samples.yaml <<~~
  > samples:
  >   A:
  >     query: region == 'North America'
  >     drop_sample: true
  >   B:
  >     target_sample: D
  >     query: country == 'Dominican Republic'
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  ERROR: Sample 'B' depends on 'D' which is not a defined sample.
  [2]

A circular DAG

  $ cat >samples.yaml <<~~
  > samples:
  >   A:
  >     target_sample: B
  >     query: region == 'North America'
  >     drop_sample: true
  >   B:
  >     target_sample: A
  >     query: country == 'Dominican Republic'
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  ERROR: Dependency cycle detected among samples.
  [2]
