Setup

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../../../bin/augur}"

Testing metadata inputs for the curate command.
Running the `passthru` subcommand since it does not do any data transformations.

Create metadata TSV file for testing.

  $ cat >$TMP/metadata.tsv <<~~
  > strain	country	date
  > sequence_A	USA	2020-10-01
  > sequence_B	USA	2020-10-02
  > sequence_C	USA	2020-10-03
  > ~~

Test TSV metadata input

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.tsv
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}

Create metadata CSV file for testing.

  $ cat >$TMP/metadata.csv <<~~
  > strain,country,date
  > sequence_A,USA,2020-10-01
  > sequence_B,USA,2020-10-02
  > sequence_C,USA,2020-10-03
  > ~~

Test CSV metadata input

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.csv
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
