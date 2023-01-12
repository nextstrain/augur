Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Import metadata.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2000-01-01
  > SEQ2	2000-02-01
  > SEQ3	2000-03-01
  > SEQ4	2000-04-01
  > ~~

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv \
  >  --output-db data.sqlite3

Import the same data into another table.

  $ ${AUGUR} db import \
  >  --metadata metadata.tsv \
  >  --metadata-table metadata_copy \
  >  --output-db data.sqlite3

The tables in the database are identical to the original data.

  $ sqlite3 data.sqlite3 ".mode tabs" ".header on" "SELECT * FROM metadata" > out1.tsv
  $ diff metadata.tsv out1.tsv

  $ sqlite3 data.sqlite3 ".mode tabs" ".header on" "SELECT * FROM metadata_copy" > out2.tsv
  $ diff metadata.tsv out2.tsv
