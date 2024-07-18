Setup

  $ source "$TESTDIR"/_setup.sh

Testing metadata inputs for the curate command.
Running the `passthru` subcommand since it does not do any data transformations.

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	country	date	authors
  > sequence_A	USA	2020-10-01	A,B,C,D,E,F,G,H,I,J,K
  > sequence_B	USA	2020-10-02	A,B,C,D,E,F,G,H,I,J,K
  > sequence_C	USA	2020-10-03	A,B,C,D,E,F,G,H,I,J,K
  > ~~

Test TSV metadata input

  $ ${AUGUR} curate passthru \
  > --metadata metadata.tsv
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Test TSV metadata input from stdin

  $ cat metadata.tsv \
  >   | ${AUGUR} curate normalize-strings \
  >     --metadata -
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Create metadata CSV file for testing.

  $ cat >metadata.csv <<~~
  > strain,country,date
  > sequence_A,USA,2020-10-01
  > sequence_B,USA,2020-10-02
  > sequence_C,USA,2020-10-03
  > ~~

Test CSV metadata input

  $ ${AUGUR} curate passthru \
  > --metadata metadata.csv
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}

Test CSV metadata input from stdin

  $ cat metadata.csv \
  >   | ${AUGUR} curate normalize-strings \
  >     --metadata -
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}

Test Excel (.xls) metadata input

  $ ${AUGUR} curate passthru \
  > --metadata "$TESTDIR/../data/metadata.xls"
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Test Excel (.xlsx) metadata input

  $ ${AUGUR} curate passthru \
  > --metadata "$TESTDIR/../data/metadata.xlsx"
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Test OpenOffice (.ods) metadata input

  $ ${AUGUR} curate passthru \
  > --metadata "$TESTDIR/../data/metadata.ods"
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Excel (.xlsx) workbook, skipped rows/cols

  $ ${AUGUR} curate passthru \
  > --metadata "$TESTDIR/../data/metadata-skipped-areas.xlsx"
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Excel (.xlsx) workbook, skipped hidden sheet

  $ ${AUGUR} curate passthru \
  > --metadata "$TESTDIR/../data/metadata-skipped-hidden-sheet.xlsx"
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "authors": "A,B,C,D,E,F,G,H,I,J,K"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "authors": "A,B,C,D,E,F,G,H,I,J,K"}

Excel (.xlsx) workbook, no valid sheets

  $ ${AUGUR} curate passthru \
  > --metadata "$TESTDIR/../data/metadata-no-valid-sheet.xlsx"
  ERROR: Excel/OpenOffice workbook '*/metadata-no-valid-sheet.xlsx' contains no visible worksheets. (glob)
  
  3 other sheets found:
    - 'Hidden' (type=worksheet, visibility=hidden)
    - 'VeryHidden' (type=worksheet, visibility=veryhidden)
    - 'Chart' (type=chartsheet, visibility=visible)
  
  [2]

Create a metadata TSV file with duplicate records

  $ cat >metadata.tsv <<~~
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
  > --metadata metadata.tsv
  ERROR: Encountered record with duplicate id 'sequence_A' in .* (re)
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  [2]

Test error_all on duplicate records.

  $ ${AUGUR} curate passthru \
  > --metadata metadata.tsv \
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
  > --metadata metadata.tsv \
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
  > --metadata metadata.tsv \
  > --duplicate-reporting silent
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03"}

Test duplicate records with a bogus id column, which is expected to fail with an error.

  $ ${AUGUR} curate passthru \
  > --metadata metadata.tsv \
  > --id-column "bogus_id"
  ERROR: The provided id column 'bogus_id' does not exist in .* (re)
  [2]
