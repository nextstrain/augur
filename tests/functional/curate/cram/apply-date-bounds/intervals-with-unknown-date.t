Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

An unknown "date" will be returned as-is with `--date-field` + min/max.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --date-field date \
  >     --lower-bound cladeRootDate \
  >     --upper-bound collectionDate
  {"record": 1, "date": "2020-01-02/2020-01-23", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

An unknown "date" using min/max produces a range with the min/max dates.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --date-field date \
  >     --lower-bound cladeRootDate \
  >     --upper-bound collectionDate
  {"record": 1, "date": "2020-01-02/2020-01-23", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

An unknown "date" using only min or max leaves the date as-is.

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --date-field date \
  >     --lower-bound cladeRootDate
  {"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}

  $ echo '{"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}' \
  >   | ${AUGUR} curate apply-date-bounds \
  >     --date-field date \
  >     --upper-bound collectionDate
  {"record": 1, "date": "?", "cladeRootDate": "2020-01-02", "collectionDate": "2020-01-23"}
