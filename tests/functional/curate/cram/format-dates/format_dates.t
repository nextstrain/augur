Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing format_dates with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2020", "collectionDate": "2020-01", "releaseDate": "2020-01","updateDate": "2020-07-18T00:00:00Z"}
  > ~~

Test output with matching expected date formats.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m" "%Y-%m-%dT%H:%M:%SZ"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01-XX", "releaseDate": "2020-01-XX", "updateDate": "2020-07-18"}

Test output with unmatched expected date formats with default `ERROR_FIRST` failure reporting.
This is expected to fail with an error, so redirecting stdout since we don't care about the output.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" 1> /dev/null
  ERROR: Unable to format date string '2020-01' in field 'collectionDate' of record 0.
  [2]

Test output with unmatched expected date formats with `ERROR_ALL` failure reporting.
This is expected to fail with an error, so redirecting stdout since we don't care about the output.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" \
  >     --failure-reporting "error_all" 1> /dev/null
  ERROR: Unable to format dates for the following (record, field, date string):
  (0, 'collectionDate', '2020-01')
  (0, 'releaseDate', '2020-01')
  [2]

Test output with unmatched expected date formats while warning on failures.
This is expected to print warnings for failures and return the masked date strings for failures.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" \
  >     --failure-reporting "warn"
  WARNING: Unable to format date string '2020-01' in field 'collectionDate' of record 0.
  WARNING: Unable to format date string '2020-01' in field 'releaseDate' of record 0.
  WARNING: Unable to format dates for the following (record, field, date string):
  (0, 'collectionDate', '2020-01')
  (0, 'releaseDate', '2020-01')
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "XXXX-XX-XX", "releaseDate": "XXXX-XX-XX", "updateDate": "2020-07-18"}

Test output with unmatched expected date formats while silencing failures.
This is expected to return the masked date strings for failures.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" \
  >     --failure-reporting "silent"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "XXXX-XX-XX", "releaseDate": "XXXX-XX-XX", "updateDate": "2020-07-18"}

Test output with unmatched expected date formats while silencing failures with `--no-mask-failure`.
This is expected to return the date strings in their original format.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" \
  >     --failure-reporting "silent" \
  >     --no-mask-failure
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01", "releaseDate": "2020-01", "updateDate": "2020-07-18"}

Test output with multiple matching expected date formats.
Date with multiple matches will be parsed according to first matching format.
The "collectionDate" and "releaseDate" will match the first "%Y-%j" format, which is a complete date.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%j" "%Y-%m" "%Y-%m-%dT%H:%M:%SZ"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01-01", "releaseDate": "2020-01-01", "updateDate": "2020-07-18"}

Test output with chained format-dates commands that parses different fields with different expected formats.
Since "collectionDate" and "releaseDate" have expected formats overlap,
we can split them into two chained commands that parses them with different expected formats to produce the desired results.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m" "%Y-%m-%dT%H:%M:%SZ" \
  >   | ${AUGUR} curate format-dates \
  >     --date-field "collectionDate" \
  >     --expected-date-formats "%Y-%j"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01-01", "releaseDate": "2020-01-XX", "updateDate": "2020-07-18"}
