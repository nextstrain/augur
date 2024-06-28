Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Overwriting an existing abbr-author-field generates a warning

`--abbr-authors-field` can be used to put the abbreviated authors into a different field

  $ echo '{"authors":"Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A.", "abbr-authors":"I EXIST"}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --abbr-authors-field="abbr-authors"
  WARNING: the 'abbr-authors' field already exists in record 0 and will be overwritten!
  {"authors": "Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A.", "abbr-authors": "Troesemeier et al."}
