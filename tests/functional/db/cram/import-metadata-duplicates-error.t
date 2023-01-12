Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Metadata with a duplicate under the ID column will error.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2000-01-01
  > SEQ1	2000-02-01
  > SEQ3	2000-03-01
  > ~~

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv \
  >  --output-db data.sqlite3
  ERROR: Duplicate found in column 'strain' of table 'metadata'.
  [2]
