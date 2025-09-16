Setup

  $ source "$TESTDIR"/_setup.sh

Create an include file.

  $ cat >include.txt <<~~
  > EcEs062_16
  > ~~

'include' can be specified as a list:

  $ cat >config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - include.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config.yaml'...
  9 strains were dropped during filtering
  	1 was added back because it was in include.txt
  	10 were dropped because of subsampling criteria
  3 strains passed all filters

'include' can be specified as a scalar:

  $ cat >config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include: include.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config.yaml'...
  9 strains were dropped during filtering
  	1 was added back because it was in include.txt
  	10 were dropped because of subsampling criteria
  3 strains passed all filters
