Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test an exact rule with duplicate raw entries applies the last matching annotation.
Verifies that we can override "default" rules by appending custom rules.

  $ cat >rules.tsv <<~~
  > r_raw/c_raw/d_raw/l_raw	r_annotated/c_annotated/d_annotated/l_annotated
  > r_raw/c_raw/d_raw/l_raw	r_annotated_custom/c_annotated_custom/d_annotated_custom/l_annotated_custom
  > ~~

It applies the last rules that matched.

  $ echo '{"region": "r_raw", "country": "c_raw", "division": "d_raw", "location": "l_raw"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  {"region": "r_annotated_custom", "country": "c_annotated_custom", "division": "d_annotated_custom", "location": "l_annotated_custom"}
