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

  $ cat >z.fasta <<~~
  > >one
  > ATCG
  > >two
  > GCTA
  > >three
  > GACT
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences X=x.fasta Y=y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta
  Reading 'X' metadata from 'x.tsv'…
  Reading 'Y' metadata from 'y.tsv'…
  Reading sequence IDs from 'x.fasta'…
  Reading sequence IDs from 'y.fasta'…
  Merging metadata and writing to 'merged.tsv'…
  Reading sequences from 'y.fasta'…
  Reading sequences from 'x.fasta'…
  Merging sequences and writing to 'merged.fasta'…
  [SeqKit] [INFO]\x1b[0m 1 duplicated records removed (esc)

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

Sequence names are optional.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv Z=z.tsv \
  >   --sequences x.fasta y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet

If sequence names are provided,

(1) they must have matching metadata inputs.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences X=x.fasta Y=y.fasta Z=z.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet
  ERROR: Named inputs must be paired.
  
  The following sequence input does not have a corresponding metadata input:
  
    'z.fasta'
  
  [2]

(2) they must exhaustively cover all metadata inputs.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv Z=z.tsv \
  >   --sequences X=x.fasta Y=y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet
  ERROR: Named inputs must be paired.
  
  The following metadata input does not have a corresponding sequence input:
  
    'z.tsv'
  
  [2]

(3) they must be in the same order as metadata inputs.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences Y=y.fasta X=x.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet
  ERROR: Named inputs must be paired with the same ordering.
  
  Order of inputs differs between named metadata ['X', 'Y'] and named sequences ['Y', 'X'].
  
  [2]

(4) IDs are cross-checked between paired input files, but only a warning is emitted.

  $ cat >>y.fasta <<~~
  > >four
  > ATCG
  > ~~

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences X=x.fasta Y=y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta
  Reading 'X' metadata from 'x.tsv'\xe2\x80\xa6 (esc)
  Reading 'Y' metadata from 'y.tsv'\xe2\x80\xa6 (esc)
  Reading sequence IDs from 'x.fasta'\xe2\x80\xa6 (esc)
  Reading sequence IDs from 'y.fasta'\xe2\x80\xa6 (esc)
  WARNING: Sequence 'four' in 'y.fasta' is missing from 'y.tsv'. Outputs may continue to be mismatched.
  Merging metadata and writing to 'merged.tsv'\xe2\x80\xa6 (esc)
  Reading sequences from 'y.fasta'\xe2\x80\xa6 (esc)
  Reading sequences from 'x.fasta'\xe2\x80\xa6 (esc)
  Merging sequences and writing to 'merged.fasta'\xe2\x80\xa6 (esc)
  [SeqKit] [INFO]\x1b[0m 1 duplicated records removed (esc)

(5) Unnamed files can still be present anywhere in the list.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv Y=y.tsv \
  >   --sequences X=x.fasta z.fasta Y=y.fasta \
  >   --output-metadata merged.tsv \
  >   --output-sequences merged.fasta \
  >   --quiet
