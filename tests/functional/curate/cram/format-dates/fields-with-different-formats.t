Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing format_dates with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2020", "collectionDate": "2020-01", "releaseDate": "2020-01","updateDate": "2020-07-18T00:00:00Z"}
  > ~~

Test output with matching expected date formats for multiple fields.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m" "%Y-%m-%dT%H:%M:%SZ"
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01-XX", "releaseDate": "2020-01-XX", "updateDate": "2020-07-18"}


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
