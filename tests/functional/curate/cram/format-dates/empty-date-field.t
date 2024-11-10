Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test empty date value, which should be returned as a fully masked date.

  $ echo '{"record": 1, "date": ""}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date"
  {"record": 1, "date": "XXXX-XX-XX"}

Test whitespace only date value, which should be returned as a fully masked date.

  $ echo '{"record": 1, "date": "  "}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date"
  {"record": 1, "date": "XXXX-XX-XX"}
