SETUP

  $ source "$TESTDIR"/_setup.sh


BASIC USAGE

Both metadata and sequences are merged with duplicate ID handling.

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

  $ cat >z.tsv <<~~
  > strain	g	c
  > one	Z1g	
  > two	Z2g	Z2c
  > three	Z3g	
  > ~~

  $ cat >x.fasta <<~~
  > >one
  > ATCG
  > >two
  > GCTA
  > ~~

  $ cat >y.fasta <<~~
  > >two
  > ATCG
  > >three
  > GCTA
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences x.fasta y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta
  Reading 'X' metadata from 'x.tsv'…
  Reading 'Y' metadata from 'y.tsv'…
  Merging metadata and writing to 'merged.tsv'…
  Validating 'x.fasta'…
  Validating 'y.fasta'…
  Merging sequences and writing to 'merged.fasta'…

  $ tsv-pretty merged.tsv
  strain  a    b    c    f    e    d
  one     X1a  X1b  X1c            
  two     X2a  X2b  Y2c  Y2f  Y2e  Y2d
  three                  Y3f  Y3e  Y3d

  $ cat merged.fasta
  >two
  ATCG
  >three
  GCTA
  >one
  ATCG

Sequence inputs are unnamed, so what seems to be a mismatch by file names alone
will be taken as-is.

  $ cat >y.fasta <<~~
  > >two
  > ATCG
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv Z=z.tsv \
  >   --sequences y.fasta x.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet

  $ cat merged.fasta
  >one
  ATCG
  >two
  GCTA
