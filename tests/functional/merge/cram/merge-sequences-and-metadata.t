SETUP

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Merge sequences and metadata

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

  $ cat >x.tsv <<~~
  > strain	a	b	c
  > one	X1a	X1b	X1c
  > two	X2a	X2b	X2c
  > ~~

  $ cat >y.tsv <<~~
  > strain	b	c	f	e	d
  > two		Y2c	Y2f	Y2e	Y2d
  > three			Y3f	Y3e	Y3d
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences X=x.fasta Y=y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet
  [INFO]\x1b[0m 1 duplicated records removed (esc)
