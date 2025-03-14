Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

An ambiguous "date" will be limited by "collectionDate".

  $ echo '{"record": 1, "date": "2020", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "2020-01-01/2020-01-23", "collectionDate": "2020-01-23"}

If "date" is an exact date, it will be left as-is.

  $ echo '{"record": 1, "date": "2020-01-10", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "2020-01-10", "collectionDate": "2020-01-23"}
