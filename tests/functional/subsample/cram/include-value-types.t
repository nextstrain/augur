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
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/include.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters

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
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/include.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters
