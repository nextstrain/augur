Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing.

  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2021", "collectionDate": "2020-01-23"}
  > {"record": 2, "date": "2022", "collectionDate": "2020-01-23"}
  > ~~

The default behavior of data error handling is to stop on the first error.

  $ cat records.ndjson \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --date-field date \
  >     --upper-bound collectionDate 1> /dev/null
  ERROR: [record 0] 'date'='2021' is later than the upper bound of 'collectionDate'='2020-01-23'
  [2]

Data errors can be batch reported all at once.

  $ cat records.ndjson \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --failure-reporting "error_all" \
  >     --date-field date \
  >     --upper-bound collectionDate 1> /dev/null
  ERROR: Unable to apply bounds. All errors:
  [record 0] 'date'='2021' is later than the upper bound of 'collectionDate'='2020-01-23'
  [record 1] 'date'='2022' is later than the upper bound of 'collectionDate'='2020-01-23'
  [2]

Data errors can emit warnings instead of a failure.

  $ cat records.ndjson \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --failure-reporting "warn" \
  >     --date-field date \
  >     --upper-bound collectionDate
  WARNING: [record 0] 'date'='2021' is later than the upper bound of 'collectionDate'='2020-01-23'
  WARNING: [record 1] 'date'='2022' is later than the upper bound of 'collectionDate'='2020-01-23'

Errors regarding the bounds themselves are not considered data errors and will stop on
the first error regardless of --failure-reporting.

  $ cat records.ndjson \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --failure-reporting "silent" \
  >     --date-field date \
  >     --upper-bound collectionDate2
  ERROR: Expected --upper-bound to be a field name or date, but got 'collectionDate2'.
  [2]
