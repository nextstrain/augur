Setup

  $ source "$TESTDIR"/_setup.sh

File paths in the config must be relative to the working directory, not the
location of the config file.

  $ mkdir -p defaults/
  $ cat >defaults/include.txt <<~~
  > EcEs062_16
  > ~~

  $ mkdir -p config/
  $ cat >config/config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - defaults/include.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  9 strains were dropped during filtering
  	1 was added back because it was in defaults/include.txt
  	10 were dropped because of subsampling criteria
  3 strains passed all filters
