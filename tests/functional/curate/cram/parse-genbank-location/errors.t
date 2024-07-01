Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Records without a `database` field result in the expected warning

  $ echo '{"geo_loc_name":"Canada:Vancouver"}' \
  >   | ${AUGUR} curate parse-genbank-location
  Record must contain `database` field to use `transform-genbank-location.`
  {"geo_loc_name": "Canada:Vancouver"}

Records with a `database` field with an unsupported value result in the expected warning

  $ echo '{"geo_loc_name":"Canada:Vancouver", "database":"database"}' \
  >   | ${AUGUR} curate parse-genbank-location
  Database value of database not supported for `transform-genbank-location`; must be "GenBank" or "RefSeq".
  {"geo_loc_name": "Canada:Vancouver", "database": "database"}

Records without a `location` field result in the expected warning

  $ echo '{"database":"GenBank"}' \
  >   | ${AUGUR} curate parse-genbank-location
  `parse-genbank-location` requires a `geo_loc_name` field; this record does not have one.
  {"database": "GenBank"}

Records without a `location` field and a custom `--location-field` result in the expected warning

  $ echo '{"database":"GenBank"}' \
  >   | ${AUGUR} curate parse-genbank-location --location-field location
  `parse-genbank-location` requires a `location` field; this record does not have one.
  {"database": "GenBank"}
