Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

We sometimes use question marks '?' as an unknown field. `augur curate apply-record-annotations`
special-cases these so that they aren't interpreted as a valid value.

This test uses real data from seasonal-flu where this bug was first observed

Rules taken directly from `augur/data/geolocation_rules.tsv` as at c8181c7

  $ cat >rules.tsv <<~~
  > Asia/Bangladesh/Jashore/	Asia/Bangladesh/Khulna/Jashore
  > Asia/Bangladesh/Jashore/*	Asia/Bangladesh/Khulna/*
  > ~~

Fist check that if the division is EMPTY then the first rule is applied, as we expect

  $ echo '{"region": "Asia", "country": "Bangladesh", "division": "Jashore", "location": ""}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "Asia", "country": "Bangladesh", "division": "Khulna", "location": "Jashore"}

but if the location is a "?", which is the case in our seasonal-flu ingest pipeline we want to use the first rule too
rather than interpreting "?" as a valid location to preserve


  $ echo '{"region": "Asia", "country": "Bangladesh", "division": "Jashore", "location": "?"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "Asia", "country": "Bangladesh", "division": "Khulna", "location": "Jashore"}


Edge case: preserve '?' in the output if we don't replace it via a rule

  $ echo '{"region": "Asia", "country": "X", "division": "?", "location": "?"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "Asia", "country": "X", "division": "?", "location": "?"}


This works for multiple fields in a work-backwards fashion

  $ cat >rules.tsv <<~~
  > North America/USA Alabama//	North America/USA/Alabama/
  > ~~

  $ echo '{"region": "North America", "country": "USA Alabama", "division": "", "location": ""}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "North America", "country": "USA", "division": "Alabama", "location": ""}

  $ echo '{"region": "North America", "country": "USA Alabama", "division": "?", "location": "?"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "North America", "country": "USA", "division": "Alabama", "location": ""}

NOTE: but division=? and location=empty string (or vice-versa) doesn't work - it's expected that the missing-value is used consistently
  $ echo '{"region": "North America", "country": "USA Alabama", "division": "?", "location": ""}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "North America", "country": "USA Alabama", "division": "?", "location": ""}

  $ echo '{"region": "North America", "country": "USA Alabama", "division": "", "location": "?"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules --geolocation-rules rules.tsv
  {"region": "North America", "country": "USA Alabama", "division": "", "location": "?"}
