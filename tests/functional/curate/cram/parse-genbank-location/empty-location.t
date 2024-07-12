Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"
  $ echo -e '{"geo_loc_name":"Canada:Vancouver"}\n{"geo_loc_name":""}' > test.ndjson

Parsing multiple records, some of which have an empty location field, works as expected.

  $  ${AUGUR} curate parse-genbank-location < test.ndjson
  {"geo_loc_name": "Canada:Vancouver", "country": "Canada", "division": "Vancouver", "location": ""}
  {"geo_loc_name": "", "country": "", "division": "", "location": ""}
