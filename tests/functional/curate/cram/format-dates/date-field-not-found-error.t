Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Providing a date field that does not exist in the record should result in an error.

  $ echo '{"record": 1, "date": "2024-01-01"}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "bad-date-field"
  ERROR: Expected date field 'bad-date-field' not found in record 0.
  [2]
