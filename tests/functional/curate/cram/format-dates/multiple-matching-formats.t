Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing format_dates with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2020", "collectionDate": "2020-01", "releaseDate": "2020-01","updateDate": "2020-07-18T00:00:00Z"}
  > ~~

Test output with multiple matching expected date formats.
Date with multiple matches will be parsed according to first matching format.
The "collectionDate" and "releaseDate" will match the first "%Y-%j" format, which is a complete date.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%j" "%Y-%m" "%Y-%m-%dT%H:%M:%SZ"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01-01", "releaseDate": "2020-01-01", "updateDate": "2020-07-18"}
