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

An unknown "date" using a combination of `--date-fields` + `--target-date-field` + only min or max leaves the date as-is.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields date \
  >     --expected-date-formats '?' \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate
  {"record": 1, "date": "XXXX-XX-XX", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --date-fields date \
  >     --expected-date-formats '?' \
  >     --target-date-field date \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "XXXX-XX-XX", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}
