Setup

  $ source "$TESTDIR"/_setup.sh

Test that subsampling config is validated against schema.

  $ cat >valid_samples.yaml <<~~
  > samples:
  >   test_sample:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 1
  >     exclude_where:
  >     - region!=South America
  > ~~

Valid config should pass schema validation.

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
  12 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t1 was added back because it was in .*sample_test_sample.* (re)
  1 strain passed all filters

Test config-root validation works for nested configs.

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

Nested config with --config-root should also pass validation.

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
  12 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t1 was added back because it was in .*sample_test_sample.* (re)
  1 strain passed all filters
