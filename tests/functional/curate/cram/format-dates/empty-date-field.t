Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test empty date value.
This currently has the unexpected behavior of returning the empty string.

  $ echo '{"record": 1, "date": ""}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date"
  {"record": 1, "date": ""}

Test whitespace only date value.
This currently raises an error.

  $ echo '{"record": 1, "date": "  "}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date"
  ERROR: Unable to format date string '  ' in field 'date' of record 0.
  [2]
