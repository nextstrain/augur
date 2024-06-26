Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Using invalid rules results in warning messages.
Valid rule matches are still applied.

  $ cat >rules.tsv <<~~
  > r_old/c_old/d_old/l_old
  > r_old/c_old/d_old/l_old	r_new/c_new/d_new
  > r_old/c_old/d_old	r_new/c_new/d_new/d_new
  > r_old/c_old/d_old/l_old	r_new/c_new/d_new/l_new
  > ~~

  $ echo '{"region": "r_old", "country": "c_old", "division": "d_old", "location": "l_old"}' \
  >   |  ${AUGUR} curate apply-geolocation-rules \
  >       --geolocation-rules rules.tsv
  WARNING: Could not decode geolocation rule 'r_old/c_old/d_old/l_old\n'. Please make sure rules are formatted as 'region/country/division/location<tab>region/country/division/location'.
  WARNING: Could not decode the annotated geolocation 'r_new/c_new/d_new'. Please make sure it is formatted as 'region/country/division/location'.
  WARNING: Could not decode the raw geolocation 'r_old/c_old/d_old'. Please make sure it is formatted as 'region/country/division/location'.
  {"region": "r_new", "country": "c_new", "division": "d_new", "location": "l_new"}
