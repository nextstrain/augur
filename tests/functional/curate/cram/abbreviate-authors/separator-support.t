Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Semi-colon separator is supported

  $ echo '{"authors":"Troesemeier,J.-H.; Musso,D.; Bluemel,J.; and Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}

Ampersand separator is supported

  $ echo '{"authors":"Troesemeier,J.-H., Musso,D., Bluemel,J. & Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}

Semi-colons and ampersand separators together are supported

  $ echo '{"authors":"Troesemeier,J.-H.; Musso,D.; Bluemel,J. & Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}
