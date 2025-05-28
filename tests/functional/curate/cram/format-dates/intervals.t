Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

An ambiguous "date" will be limited by a lower and upper bound based on other fields.

  $ echo '{"record": 1, "date": "2020", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "2020-01-02/2020-01-23", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

It can be limited by a lower bound only.

  $ echo '{"record": 1, "date": "2020", "cladeRootDate": "2020-01-02"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate
  {"record": 1, "date": "2020-01-02/2020-12-31", "cladeRootDate": "2020-01-02"}

It can be limited by an upper bound only.

  $ echo '{"record": 1, "date": "2020", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "2020-01-01/2020-01-23", "collectionDate": "2020-01-23"}

An error is shown if it cannot fall within bounds.

  $ echo '{"record": 1, "date": "2019", "cladeRootDate": "2020-01-02"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate
  ERROR: date='2019' is earlier than the lower bound of cladeRootDate='2020-01-02'
  [2]

  $ echo '{"record": 1, "date": "2021", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-max collectionDate
  ERROR: date='2021' is later than the upper bound of collectionDate='2020-01-23'
  [2]

Exact dates are kept as-is.

  $ echo '{"record": 1, "date": "2020-01-10", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate format-dates \
  >     --target-date-field date \
  >     --target-date-field-min cladeRootDate \
  >     --target-date-field-max collectionDate
  {"record": 1, "date": "2020-01-10", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}
