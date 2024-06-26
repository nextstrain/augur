Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that closer wildcard matches (with respect to hierarchy) are applied before
more generic wildcard matches.

  $ cat >rules.tsv <<~~
  > a/b/c/*	1/2/3/*
  > a/b/*/*	1/2/*/*
  > a/*/*/*	1/*/*/*
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "x"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "2", "division": "3", "location": "x"}

  $ echo '{"region": "a", "country": "b", "division": "x", "location": "x"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "2", "division": "x", "location": "x"}

  $ echo '{"region": "a", "country": "x", "division": "x", "location": "x"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "x", "division": "x", "location": "x"}
