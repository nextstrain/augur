Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

If the authors list already ends in `et al.` don't add another.

  $ echo '{"authors":"Troesemeier et al."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}
