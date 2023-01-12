Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Metadata with a row that has an extra column does not error.
A missing column also does not error.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2000-01-01
  > SEQ2	2000-02-01	extra
  > SEQ3
  > ~~

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv \
  >  --output-db data.sqlite3

The extra column is ignored, and missing columns are empty.

  $ sqlite3 data.sqlite3 ".mode tabs" ".header on" "SELECT * FROM metadata"
  strain\tdate (esc)
  SEQ1\t2000-01-01 (esc)
  SEQ2\t2000-02-01 (esc)
  SEQ3\t (esc)
