Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Import sequences.

  $ cat >sequences.fasta <<~~
  > > SEQ1
  > ATAT
  > > SEQ2
  > CGGC
  > > SEQ3
  > ACGT
  > ~~

  $ ${AUGUR} db import \
  >  --sequences sequences.fasta \
  >  --output-db data.sqlite3

Import the same data into another table.

  $ ${AUGUR} db import \
  >  --sequences sequences.fasta \
  >  --sequences-table sequences_copy \
  >  --output-db data.sqlite3

The tables in the database are identical to the original data.

  $ sqlite3 data.sqlite3 ".mode csv" ".header on" "SELECT * FROM sequences"
  strain,sequence\r (esc)
  SEQ1,ATAT\r (esc)
  SEQ2,CGGC\r (esc)
  SEQ3,ACGT\r (esc)

  $ sqlite3 data.sqlite3 ".mode csv" ".header on" "SELECT * FROM sequences_copy"
  strain,sequence\r (esc)
  SEQ1,ATAT\r (esc)
  SEQ2,CGGC\r (esc)
  SEQ3,ACGT\r (esc)
