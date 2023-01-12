Setup

  $ export AUGUR="$TESTDIR/../../../../bin/augur"

Import metadata and sequences.

  $ cat >in.tsv <<~~
  > strain	date
  > SEQ1	2000-01-01
  > SEQ2	2000-02-01
  > SEQ3	2000-03-01
  > SEQ4	2000-04-01
  > ~~

  $ cat >in.fasta <<~~
  > > SEQ1
  > ATAT
  > > SEQ2
  > CGGC
  > > SEQ3
  > ACGT
  > ~~

  $ ${AUGUR} db import \
  >  --metadata in.tsv \
  >  --sequences in.fasta \
  >  --output-db data.sqlite3

Query and create new tables in the database.

  $ sqlite3 data.sqlite3 "
  >   CREATE TABLE filtered_metadata AS
  >   SELECT * FROM metadata
  >   WHERE strain in ('SEQ1', 'SEQ3')
  > "
  $ sqlite3 data.sqlite3 "
  >   CREATE TABLE filtered_sequences AS
  >   SELECT * FROM sequences
  >   WHERE strain in ('SEQ1', 'SEQ3')
  > "

Export the new tables.

  $ ${AUGUR} db export \
  >  --db data.sqlite3 \
  >  --metadata-table filtered_metadata \
  >  --sequences-table filtered_sequences \
  >  --output-metadata out.tsv \
  >  --output-sequences out.fasta
  $ cat out.tsv
  strain\tdate (esc)
  SEQ1\t2000-01-01 (esc)
  SEQ3\t2000-03-01 (esc)
  $ cat out.fasta
  >SEQ1
  ATAT
  >SEQ3
  ACGT
