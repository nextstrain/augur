Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

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
