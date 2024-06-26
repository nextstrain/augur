Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Create a file with cyclic rules.

  $ cat >rules.tsv <<~~
  > r_old/c_old/d_old/l_old	r_new/c_new/d_new/l_new
  > r_new/c_new/d_new/l_new	r_old/c_old/d_old/l_old
  > ~~

Attempting to use the rules with no matches results in no changes.

  $ echo '{"region": "r_old", "country": "c_old", "division": "d_old", "location": "l_something_else"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "r_old", "country": "c_old", "division": "d_old", "location": "l_something_else"}

Attempting to use the rules with a match results in an error.

  $ echo '{"region": "r_old", "country": "c_old", "division": "d_old", "location": "l_old"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  ERROR: More than 1000 geolocation rules applied on the same entry ['r_old', 'c_old', 'd_old', 'l_old'].
  [2]
