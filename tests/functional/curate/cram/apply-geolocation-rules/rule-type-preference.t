Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that matching rules are used over general rules.

  $ cat >rules.tsv <<~~
  > a/b/c/d	1/2/3/4
  > a/*/*/*	2/*/*/*
  > ~~

  $ echo '{"region": "a", "country": "b", "division": "c", "location": "d"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "1", "country": "2", "division": "3", "location": "4"}
