Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing format_dates with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2020", "collectionDate": "2020-01", "releaseDate": "2020-01","updateDate": "2020-07-18T00:00:00Z"}
  > ~~

Test output with unmatched expected date formats while silencing failures with `--no-mask-failure`.
This is expected to return the date strings in their original format.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" "collectionDate" "releaseDate" "updateDate" \
  >     --expected-date-formats "%Y" "%Y-%m-%dT%H:%M:%SZ" \
  >     --failure-reporting "silent" \
  >     --no-mask-failure
  {"record": 1, "date": "2020-XX-XX", "collectionDate": "2020-01", "releaseDate": "2020-01", "updateDate": "2020-07-18"}
