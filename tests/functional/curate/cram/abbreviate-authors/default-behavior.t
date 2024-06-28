Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Running the command with no arguments produces the expected output

  $ echo '{"authors":"Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}

`--authors-field` can be used to set an alternative field name

  $ echo '{"author-list":"Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --authors-field="author-list"
  {"author-list": "Troesemeier et al."}

`--default-value` can be used to provide a default for an empty field

  $ echo '{"authors":""}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --default="??"
  {"authors": "??"}

`--abbr-authors-field` can be used to put the abbreviated authors into a different field

  $ echo '{"authors":"Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --abbr-authors-field="abbr-authors"
  {"authors": "Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A.", "abbr-authors": "Troesemeier et al."}

`--authors-field` and `--default-value` work together

  $ echo '{"author-list":""}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --authors-field="author-list" \
  >       --default-value="???"
  {"author-list": "???"}

`--authors-field` and `--abbr-authors-field` work together

  $ echo '{"author-list":"Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --authors-field="author-list" \
  >       --abbr-authors-field="abbr-authors"
  {"author-list": "Troesemeier,J.-H., Musso,D., Bluemel,J. and Baylis,S.A.", "abbr-authors": "Troesemeier et al."}

`--default-value` and `--abbr-authors-field` work together

  $ echo '{"authors":""}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --default-value="?" \
  >       --abbr-authors-field="abbr-authors"
  {"authors": "", "abbr-authors": "?"}

All three options work together

  $ echo '{"author-list":""}' \
  >   | ${AUGUR} curate abbreviate-authors \
  >       --authors-field="author-list" \
  >       --abbr-authors-field="abbr-authors" \
  >       --default-value="?!"
  {"author-list": "", "abbr-authors": "?!"}

Running the command with no arguments and multiple records produces the expected output

  $ echo '{"authors":"Troesemeier,J.-H. & Musso,D."}
  > {"authors":"Bluemel,J. and Baylis,S.A."}' \
  >   | ${AUGUR} curate abbreviate-authors
  {"authors": "Troesemeier et al."}
  {"authors": "Bluemel et al."}
