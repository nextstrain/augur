SETUP

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Merge sequences without metadata.

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
  >   --sequences x.fasta y.fasta \
  >   --output-sequences merged.fasta
  Validating 'x.fasta'…
  Validating 'y.fasta'…
  Merging sequences and writing to 'merged.fasta'…

  $ cat merged.fasta
  >seq3
  ATCG
  >seq4
  GCTA
  >seq1
  ATCG
  >seq2
  GCTA

Sequence files with missing newlines at the end are supported.

  $ truncate -s -1 x.fasta

  $ truncate -s -1 y.fasta

  $ ${AUGUR} merge \
  >   --sequences x.fasta y.fasta \
  >   --output-sequences merged.fasta
  Validating 'x.fasta'…
  Validating 'y.fasta'…
  Merging sequences and writing to 'merged.fasta'…

  $ cat merged.fasta
  >seq3
  ATCG
  >seq4
  GCTA
  >seq1
  ATCG
  >seq2
  GCTA
