Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Metadata files must have at least a header row, so fail if it is empty.

  $ cat >empty.txt <<~~
  > ~~

  $ ${AUGUR} db import \
  >  --metadata empty.txt \
  >  --output-db data.sqlite3
  ERROR: Expected a header row in 'empty.txt' but it is empty.
  [2]

An empty file is a valid FASTA file, so do not error in this case.

  $ ${AUGUR} db import \
  >  --sequences empty.txt \
  >  --output-db data.sqlite3
