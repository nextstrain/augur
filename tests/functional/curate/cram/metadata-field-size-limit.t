Setup

  $ source "$TESTDIR"/_setup.sh

Testing the metadata field size limit for the curate command.
Running the `passthru` subcommand since it does not do any data transformations.

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	authors
  > sequence_A	A,B,C,D,E,F,G,H,I,J,K
  > ~~

Test field larger than the limit

  $ AUGUR_FIELD_SIZE_LIMIT=10 ${AUGUR} curate passthru \
  > --metadata metadata.tsv
  ERROR: Unable to read metadata file 'metadata.tsv': a field on line 2 exceeds the field size limit of 10 characters set by AUGUR_FIELD_SIZE_LIMIT.
  
  This is often caused by a missing closing quote. If the data is correct, rerun with AUGUR_FIELD_SIZE_LIMIT increased above 10.
  [2]

Test increased limit

  $ AUGUR_FIELD_SIZE_LIMIT=21 ${AUGUR} curate passthru \
  > --metadata metadata.tsv
  {"strain": "sequence_A", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Test invalid limit

  $ AUGUR_FIELD_SIZE_LIMIT=invalid ${AUGUR} curate passthru \
  > --metadata metadata.tsv
  ERROR: AUGUR_FIELD_SIZE_LIMIT must be a positive integer. Received 'invalid'.
  [2]
