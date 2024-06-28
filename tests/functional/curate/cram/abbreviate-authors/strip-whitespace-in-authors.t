Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Whitespace in author list gets stripped out

  $ echo '{"authors":"Troesemeier,J.-H.,    Musso,D., Bluemel,J. and     Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}
