Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

An unknown "date" will be returned as-is with `--target-date-field` + min/max.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

An unknown "date" using a combination of `--date-fields` + `--target-date-field` + min/max produces
a range with the min/max dates.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields date \
  >     --expected-date-formats '?' \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "2020-01-02/2020-01-23", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

An unknown "date" using a combination of `--date-fields` + `--target-date-field` + only min produces
a range with the max date set to today. (Note this test will fail when run on a different date.)

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields date \
  >     --expected-date-formats '?' \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate
  {"record": 1, "date": "2020-01-02/2025-04-02", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

An unknown "date" using a combination of `--date-fields` + `--target-date-field` + only max produces a range with
the min date set to `0001-01-01`.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields date \
  >     --expected-date-formats '?' \
  >     --target-date-field date \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "0001-01-01/2020-01-23", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}
