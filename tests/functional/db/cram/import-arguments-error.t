Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Missing input data.

  $ ${AUGUR} db import \
  >  --output-db data.sqlite3 2> /dev/null
  [2]

Missing output destination.

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv 2> /dev/null
  [2]
