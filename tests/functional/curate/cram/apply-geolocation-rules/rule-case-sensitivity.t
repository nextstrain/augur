Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test a rule with different casing for the annotations.

  $ cat >rules.tsv <<~~
  > North America/USA/CA/*	North America/USA/California/*
  > ~~

Rule matching is case-sensitive and the output matches the casing of the annotations.

  $ echo '{"region": "North America", "country": "USA", "division": "CA", "location": "Los Angeles"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "North America", "country": "USA", "division": "California", "location": "Los Angeles"}

Rule matching is case-sensitive, so raw values with mismatched casing do not get changed.

  $ echo '{"region": "North America", "country": "USA", "division": "Ca", "location": "Los Angeles"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "North America", "country": "USA", "division": "Ca", "location": "Los Angeles"}
