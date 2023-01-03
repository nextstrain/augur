Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Ensure the sqlite_ table name prefix error is handled, since it is forbidden¹.

¹ https://www.sqlite.org/lang_createtable.html

  $ cat >metadata.tsv <<~~
  > strain	col
  > SEQ1	A
  > SEQ2	B
  > SEQ3	C
  > ~~

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv \
  >  --metadata-table sqlite_table \
  >  --output-db data.sqlite3
  ERROR: Failed to create table: object name reserved for internal use: sqlite_table
  [2]
