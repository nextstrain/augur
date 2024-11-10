Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing format_dates with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2020", "collectionDate": "2020-01", "releaseDate": "2020-01","updateDate": "2020-07-18T00:00:00Z"}
  > ~~

Test output with unmatched expected date formats with default `ERROR_FIRST` failure reporting.
This is expected to fail with an error, so redirecting stdout since we don't care about the output.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" 1> /dev/null
  ERROR: Unable to format date string '2020-01' in field 'collectionDate' of record 0.
  Current expected date formats are ['%Y-%m-%d', '%Y-%m-XX', '%Y-XX-XX', 'XXXX-XX-XX', '%Y', '%Y-%m-%dT%H:%M:%SZ']. This can be updated with --expected-date-formats.
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
  Current expected date formats are ['%Y-%m-%d', '%Y-%m-XX', '%Y-XX-XX', 'XXXX-XX-XX', '%Y', '%Y-%m-%dT%H:%M:%SZ']. This can be updated with --expected-date-formats.
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
  Current expected date formats are ['%Y-%m-%d', '%Y-%m-XX', '%Y-XX-XX', 'XXXX-XX-XX', '%Y', '%Y-%m-%dT%H:%M:%SZ']. This can be updated with --expected-date-formats.
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "XXXX-XX-XX", "releaseDate": "XXXX-XX-XX", "updateDate": "2020-07-18"}

Test output with unmatched expected date formats while silencing failures.
This is expected to return the masked date strings for failures.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" \
  >     --failure-reporting "silent"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "XXXX-XX-XX", "releaseDate": "XXXX-XX-XX", "updateDate": "2020-07-18"}
