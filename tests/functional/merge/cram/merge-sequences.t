SETUP

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Merge sequences without metadata

  $ cat >x.fasta <<~~
  > >seq1
  > ATCG
  > >seq2
  > GCTA
  > >seq3
  > TCGA
  > ~~

  $ cat >y.fasta <<~~
  > >seq3
  > ATCG
  > >seq4
  > GCTA
  > ~~

  $ ${AUGUR} merge \
  >   --sequences x=x.fasta y=y.fasta \
  >   --output-sequences - > merged.fasta
  [INFO]\x1b[0m 1 duplicated records removed (esc)

  $ cat merged.fasta
  >seq3
  ATCG
  >seq4
  GCTA
  >seq1
  ATCG
  >seq2
  GCTA

FIXME: Merge both sequences and metadata
