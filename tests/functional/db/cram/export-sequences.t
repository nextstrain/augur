Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Create an empty sequence table.

  $ cat >in.tsv <<~~
  > strain	sequence
  > ~~

  $ ${AUGUR} db import \
  >  --metadata in.tsv \
  >  --metadata-table empty_sequences \
  >  --output-db data.sqlite3

Exporting the empty sequence table results in an empty FASTA file.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --sequences-table empty_sequences \
  >  --output-sequences out.fasta
  $ cat out.fasta

Create an invalid sequence table.

  $ cat >in.tsv <<~~
  > strain	col2
  > SEQ1	ATAT
  > SEQ2	CGGC
  > SEQ3	ACGT
  > ~~

  $ ${AUGUR} db import \
  >  --metadata in.tsv \
  >  --metadata-table table_with_wrong_column_name \
  >  --output-db data.sqlite3

Attempting to export a table as FASTA without the right columns shows an error.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --sequences-table table_with_wrong_column_name \
  >  --output-sequences out.fasta
  ERROR: Sequence column 'sequence' missing from table 'table_with_wrong_column_name'.
  [2]
