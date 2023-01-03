Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Ensure these column names work after being sanitized:

- FILTER: this is a SQLite keyword¹.
- sqlite_col: the sqlite_ prefix is forbidden for table names².
- /col: a forward slash is know to cause problems³.

¹ https://www.sqlite.org/lang_keywords.html
² https://www.sqlite.org/lang_createtable.html
³ https://stackoverflow.com/a/3373549

  $ cat >metadata.tsv <<~~
  > strain	FILTER	sqlite_col	/col
  > SEQ1	A	A	A
  > SEQ2	B	B	B
  > SEQ3	C	C	C
  > ~~

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv \
  >  --output-db data.sqlite3
