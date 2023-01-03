Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Import metadata.

  $ cat >in.tsv <<~~
  > strain	date
  > ~~

  $ ${AUGUR} db import \
  >  --metadata in.tsv \
  >  --metadata-table empty_metadata \
  >  --output-db data.sqlite3

Export a table without any rows. The contents should be identical to the file
originally imported.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --metadata-table empty_metadata \
  >  --output-metadata out.tsv
  $ diff in.tsv out.tsv
