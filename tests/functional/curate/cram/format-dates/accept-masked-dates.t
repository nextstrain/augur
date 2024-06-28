Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create NDJSON file for testing format_dates with different forms

  $ cat >records.ndjson <<~~
  > {"record": 1, "masked-date": "2020-XX-XX", "date": "2020"}
  > ~~

Test input with masked dates works as expected if the masked format is provided.

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "masked-date" \
  >     --expected-date-formats "%Y-XX-XX"
  {"record": 1, "masked-date": "2020-XX-XX", "date": "2020"}

Test mutliple passes through format-dates where the first pass masks dates.
This is expected to fail without the masked format

  $ cat records.ndjson \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" \
  >     --expected-date-formats "%Y" \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields "date" \
  >     --expected-date-formats "%Y-%m-%d"
  ERROR: Unable to format date string '2020-XX-XX' in field 'date' of record 0.
  [2]
