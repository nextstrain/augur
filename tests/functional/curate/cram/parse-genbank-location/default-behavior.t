Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Running the command with no arguments produces the expected output

  $ echo '{"geo_loc_name":"Canada:Vancouver", "database":"GenBank"}' \
  >   | ${AUGUR} curate parse-genbank-location
  {"geo_loc_name": "Canada:Vancouver", "database": "GenBank", "country": "Canada", "division": "Vancouver", "location": ""}

`--location-field` can be used to specify a different field name

  $ echo '{"location":"Canada:Vancouver", "database":"GenBank"}' \
  >   | ${AUGUR} curate parse-genbank-location --location-field location
  {"location": "", "database": "GenBank", "country": "Canada", "division": "Vancouver"}

`RefSeq` works as a database value

  $ echo '{"geo_loc_name":"Canada:Vancouver", "database":"RefSeq"}' \
  >   | ${AUGUR} curate parse-genbank-location
  {"geo_loc_name": "Canada:Vancouver", "database": "RefSeq", "country": "Canada", "division": "Vancouver", "location": ""}

Record with only a `country` value results in expected output

  $ echo '{"geo_loc_name":"Canada", "database":"GenBank"}' \
  >   | ${AUGUR} curate parse-genbank-location
  {"geo_loc_name": "Canada", "database": "GenBank", "country": "Canada", "division": "", "location": ""}

Record with `country`, `region`, and `location` values results in expected output

  $ echo '{"geo_loc_name":"France:Cote d'\''Azur, Antibes", "database":"GenBank"}' \
  >   | ${AUGUR} curate parse-genbank-location
  {"geo_loc_name": "France:Cote d'Azur, Antibes", "database": "GenBank", "country": "France", "division": "Cote d'Azur", "location": "Antibes"}
