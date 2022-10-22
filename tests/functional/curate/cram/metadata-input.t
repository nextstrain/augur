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

Create a metadata TSV file with duplicate records

  $ cat >$TMP/metadata.tsv <<~~
  > strain	country	date
  > sequence_A	USA	2020-10-01
  > sequence_B	USA	2020-10-02
  > sequence_C	USA	2020-10-03
  > sequence_A	USA	2020-10-01
  > sequence_B	USA	2020-10-02
  > sequence_C	USA	2020-10-03
  > ~~

Test default options for duplicate records, which is expected for exit with an error on the first duplicate
There will still be output due to the nature of the chained generators in augur curate.

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.tsv
  ERROR: Encountered record with duplicate id 'sequence_A' in .* (re)
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  [2]

Test error_all on duplicate records.

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.tsv \
  > --duplicate-reporting error_all
  ERROR: The following records are duplicated in .* (re)
  'sequence_A'
  'sequence_B'
  'sequence_C'
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  [2]

Test warning on duplicate records.

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.tsv \
  > --duplicate-reporting warn
  WARNING: Encountered record with duplicate id 'sequence_A' in .* (re)
  WARNING: Encountered record with duplicate id 'sequence_B' in .* (re)
  WARNING: Encountered record with duplicate id 'sequence_C' in .* (re)
  WARNING: The following records are duplicated in .* (re)
  'sequence_A'
  'sequence_B'
  'sequence_C'
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}

Test silent on duplicate records.

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.tsv \
  > --duplicate-reporting silent
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}

Test duplicate records with a bogus id column, which is expected to fail with an error.

  $ ${AUGUR} curate passthru \
  > --metadata $TMP/metadata.tsv \
  > --id-column "bogus_id"
  ERROR: The provided id column 'bogus_id' does not exist in .* (re)
  [2]
