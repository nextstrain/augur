Setup

  $ source "$TESTDIR"/_setup.sh

Valid config passes schema validation and does not result in an error.

  $ cat >valid_samples.yaml <<~~
  > samples:
  >   test_sample:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 1
  >     exclude_where:
  >     - region!=South America
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config valid_samples.yaml \
  >   --output-metadata valid_output.tsv \
  >   --output-sequences valid_output.fasta \
  >   --subsample-seed 0
  Validating schema of 'valid_samples.yaml'...
  Sampling at 1 per group.
  11 strains were dropped during filtering
  	6 were dropped because of 'region!=South America'
  	5 were dropped because of subsampling criteria
  1 strain passed all filters

Nested config is also valid with --config-root.

  $ cat >nested_valid.yaml <<~~
  > other_tool:
  >   setting: value
  > subsample:
  >   samples:
  >     test_sample:
  >       group_by:
  >       - region
  >       subsample_max_sequences: 1
  >       exclude_where:
  >       - region!=South America
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config nested_valid.yaml \
  >   --config-root subsample \
  >   --output-metadata nested_valid_output.tsv \
  >   --output-sequences nested_valid_output.fasta \
  >   --subsample-seed 0
  Validating schema of 'nested_valid.yaml'...
  Sampling at 1 per group.
  11 strains were dropped during filtering
  	6 were dropped because of 'region!=South America'
  	5 were dropped because of subsampling criteria
  1 strain passed all filters
