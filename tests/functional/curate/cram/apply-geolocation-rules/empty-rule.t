Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that an empty field can be loaded in the raw and annotated columns.

  $ cat >rules.tsv <<~~
  > r_raw/c_raw//	r_annotated/c_annotated//
  > ~~

  $ echo '{"region": "r_raw", "country": "c_raw", "division": "", "location": ""}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "r_annotated", "country": "c_annotated", "division": "", "location": ""}

This is effectively an exact rule, so it will not apply to any record that does
not match the empty string.

  $ echo '{"region": "r_raw", "country": "c_raw", "division": "d_raw", "location": ""}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "r_raw", "country": "c_raw", "division": "d_raw", "location": ""}
