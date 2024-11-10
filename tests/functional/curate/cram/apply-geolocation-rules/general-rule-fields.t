Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that general rules work in all fields.

Field 1

  $ cat >rules.tsv <<~~
  > */b/c/d	*/2/3/4
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "a", "country": "2", "division": "3", "location": "4"}

Field 2

  $ cat >rules.tsv <<~~
  > a/*/c/d	1/*/3/4
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "b", "division": "3", "location": "4"}

Field 3

  $ cat >rules.tsv <<~~
  > a/b/*/d	1/2/*/4
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "2", "division": "c", "location": "4"}

Field 4

  $ cat >rules.tsv <<~~
  > a/b/c/*	1/2/3/*
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "2", "division": "3", "location": "d"}

Fields 2,3

  $ cat >rules.tsv <<~~
  > a/*/*/d	1/*/*/4
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "b", "division": "c", "location": "4"}

Fields 1,2,3

  $ cat >rules.tsv <<~~
  > */*/*/d	*/*/*/4
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "a", "country": "b", "division": "c", "location": "4"}
