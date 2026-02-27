Setup

  $ source "$TESTDIR"/_setup.sh

File at path relative to config file location is found

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
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/defaults/include.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters

File at path relative to current working directory is found

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
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/include.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters

File in --search-path is found

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
  >   --search-path custom_dir \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/custom_dir/include_custom.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters

Custom search path is searched before default paths (the file exists in both cwd and search_path_dir)

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
  >   --search-path search_path_dir \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/include.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters

Verify the file from search_path_dir was used (contains SG_018, not EcEs062_16 from cwd)

  $ grep "SG_018" output_metadata.tsv | cut -f1
  SG_018

Error when file not found in any search path

  $ mkdir -p config/
  $ cat >config/config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - nonexistent.txt
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --output-metadata output_metadata.tsv
  Validating schema of 'config/config.yaml'...
  ERROR: File 'nonexistent.txt' not resolvable from any of the following paths:
  
    config
    .*/include-file.t (re)
  [2]

--search-paths overrides AUGUR_SEARCH_PATHS

  $ mkdir -p env_dir/
  $ cat >env_dir/include.txt <<~~
  > EcEs062_16
  > ~~

  $ mkdir -p cli_dir/
  $ cat >cli_dir/include.txt <<~~
  > SG_018
  > ~~

  $ cat >config/config.yaml <<~~
  > samples:
  >   test:
  >     max_sequences: 2
  >     include:
  >     - include.txt
  > ~~

  $ AUGUR_SEARCH_PATHS=env_dir ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config/config.yaml \
  >   --search-paths cli_dir \
  >   --output-metadata output_metadata.tsv \
  >   --seed 0
  Validating schema of 'config/config.yaml'...
  WARNING: Both the command line argument --search-paths
  and the environment variable AUGUR_SEARCH_PATHS are set.
  Only the command line argument will be used.
  [test] running as the only filter call necessary
  [test] 9 strains were dropped during filtering
  \[test\] \\t1 was added back because it was in .*/cli_dir/include.txt.* (re)
  [test] 	10 were dropped because of subsampling criteria
  [test] 3 strains passed all filters

Verify the file from cli_dir was used (contains SG_018, not EcEs062_16 from env_dir)

  $ grep "SG_018" output_metadata.tsv | cut -f1
  SG_018
