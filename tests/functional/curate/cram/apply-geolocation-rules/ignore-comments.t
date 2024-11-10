Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that comments are ignored.

  $ cat >rules.tsv <<~~
  > # This is a comment.
  > r_raw/c_raw/d_raw/l_raw	r_annotated/c_annotated/d_annotated/l_annotated # This is also a comment.
  > ~~

  $ echo '{"region": "r_raw", "country": "c_raw", "division": "d_raw", "location": "l_raw"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "r_annotated", "country": "c_annotated", "division": "d_annotated", "location": "l_annotated"}
