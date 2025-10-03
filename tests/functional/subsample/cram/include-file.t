Setup

  $ source "$TESTDIR"/_setup.sh

Test 1: File at path relative to config file location is found

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
  >     - ../defaults/include.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  9 strains were dropped during filtering
  \\t1 was added back because it was in .*/defaults/include.txt.* (re)
  	10 were dropped because of subsampling criteria
  3 strains passed all filters

Test 2: File at path relative to current working directory is found

  $ cat >include_cwd.txt <<~~
  > EcEs062_16
  > ~~

  $ cat >config/config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - include_cwd.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  9 strains were dropped during filtering
  \\t1 was added back because it was in .*/include_cwd.txt.* (re)
  	10 were dropped because of subsampling criteria
  3 strains passed all filters

Test 3: File in --config-filepath-location is found

  $ mkdir -p custom_dir/
  $ cat >custom_dir/include_custom.txt <<~~
  > EcEs062_16
  > ~~

  $ cat >config/config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - include_custom.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --config-filepath-location custom_dir \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  9 strains were dropped during filtering
  \\t1 was added back because it was in .*/custom_dir/include_custom.txt.* (re)
  	10 were dropped because of subsampling criteria
  3 strains passed all filters

Test 4: --config-filepath-location overrides default behavior (the file exists in both cwd and search_path_dir)

  $ mkdir -p search_path_dir/

  $ cat >config/include.txt <<~~
  > EcEs062_16
  > ~~

  $ cat >search_path_dir/include.txt <<~~
  > SG_018
  > ~~

  $ cat >config/config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - include.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --config-filepath-location search_path_dir \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  9 strains were dropped during filtering
  \\t1 was added back because it was in .*/include.txt.* (re)
  	10 were dropped because of subsampling criteria
  3 strains passed all filters

Verify the file from search_path_dir was used (contains SG_018, not EcEs062_16 from cwd)

  $ grep "SG_018" output_metadata.tsv | cut -f1
  SG_018

Test 5: Error when file not found in any search path

  $ mkdir -p config/
  $ cat >config/config5.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - nonexistent.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config5.yaml \
  >   --output-metadata output_metadata.tsv 2>&1 | grep "not found"
  ERROR: File 'nonexistent.txt' not found in any of the search paths.* (re)
  [2]
