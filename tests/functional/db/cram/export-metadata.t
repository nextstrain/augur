Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Import metadata.

  $ cat >in.tsv <<~~
  > strain	date
  > SEQ1	2000-01-01
  > SEQ2	2000-02-01
  > SEQ3	2000-03-01
  > SEQ4	2000-04-01
  > ~~

  $ ${AUGUR} db import \
  >  --metadata in.tsv \
  >  --output-db data.sqlite3

Exported metadata file is identical to the original imported file.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --output-metadata out.tsv
  $ diff in.tsv out.tsv
