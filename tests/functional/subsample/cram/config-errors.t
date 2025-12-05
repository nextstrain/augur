Setup

  $ source "$TESTDIR"/_setup.sh

Test handling of YAML syntax errors.

  $ cat >config.yaml <<~~
  > samples:
  >   test:
  >     invalid: [
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata output.tsv
  ERROR: The configuration file 'config.yaml' is not valid YAML.
  while parsing a flow node
  expected the node content, but found '<stream end>'
    in "config.yaml", line 4, column 1
  [2]

Test handling of schema validation errors with invalid config.

  $ cat >config.yaml <<~~
  > subsample:
  >   samples:
  >     max_sequences: 5
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --config config.yaml \
  >   --output-metadata output.tsv
  Validating schema of 'config.yaml'...
    top level failed: Unexpected property 'subsample'
    top level failed: Missing required property 'samples'
  ERROR: Validation of 'config.yaml' failed.
  [2]
