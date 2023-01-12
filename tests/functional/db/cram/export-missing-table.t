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

Attempting to export a missing metadata table shows an error.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --metadata-table missing_table \
  >  --output-metadata out.tsv
  ERROR: Table 'missing_table' does not exist.
  [2]

Attempting to export a missing sequences table shows an error.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --sequences-table missing_table \
  >  --output-sequences out.fasta
  ERROR: Table 'missing_table' does not exist.
  [2]
