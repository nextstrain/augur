SETUP

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Duplicates within the same file are not allowed.

  $ cat >x.fasta <<~~
  > >seq1
  > ATCG
  > >seq2
  > GCTA
  > >seq3
  > TCGA
  > >seq3
  > TCGATC
  > >seq3
  > TCGATCGA
  > ~~

  $ cat >y.fasta <<~~
  > >seq3
  > ATCG
  > >seq4
  > GCTA
  > ~~

  $ ${AUGUR} merge \
  >   --sequences x.fasta y.fasta \
  >   --output-sequences - > merged.fasta
  Validating 'x.fasta'…
  ERROR: File x.fasta has duplicate IDs: 'seq3'
  [2]

  $ cat merged.fasta

Duplicate checking can be skipped, however it will lead to dropping of all but
the first entry, which is at odds with the approach used for the rest of
augur merge. This behavior is a reflection of the seqkit rmdup internals.

  $ ${AUGUR} merge \
  >   --sequences x.fasta y.fasta \
  >   --output-sequences merged.fasta \
  >   --skip-input-sequences-validation \
  >   --quiet

  $ cat merged.fasta
  >seq3
  ATCG
  >seq4
  GCTA
  >seq1
  ATCG
  >seq2
  GCTA
