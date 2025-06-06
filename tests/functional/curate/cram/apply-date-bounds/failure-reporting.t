Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing.


  $ cat >records.ndjson <<~~
  > {"record": 1, "date": "2021", "collectionDate": "2020-01-23"}
  > {"record": 2, "date": "2022", "collectionDate": "2020-01-23"}
  > ~~

Test output with date out of range with "error_all" failure reporting and upper bound.

  $ cat records.ndjson \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --failure-reporting "error_all" \
  >     --date-field date \
  >     --upper-bound collectionDate 1> /dev/null
  ERROR: Unable to apply bounds for the following (record, date, lower bound, upper bound):
  (0, '2021', None, 2020.061475409836)
  (1, '2022', None, 2020.061475409836)
  The upper bound is 'collectionDate'. This can be updated with --upper-bound.
  [2]
