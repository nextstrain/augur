SETUP

  $ source "$TESTDIR"/_setup.sh

BASIC USAGE

Sequence inputs are merged with duplicate ID handling.

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
  >   --sequences x.fasta y.fasta \
  >   --output-sequences - > merged.fasta
  Reading sequence IDs from 'x.fasta'…
  Reading sequence IDs from 'y.fasta'…
  Reading sequences from 'y.fasta'…
  Reading sequences from 'x.fasta'…
  Merging sequences and writing to '-'…
  [SeqKit] [INFO]\x1b[0m 1 duplicated records removed (esc)

seq3 is in both x and y. It is taken from the latter.

  $ cat merged.fasta
  >two
  ATCG
  >three
  GCTA
  >one
  ATCG

Duplicates are not allowed within individual sequence inputs.

  $ cat >>y.fasta <<~~
  > >three
  > ATCG
  > ~~

  $ ${AUGUR} merge \
  >   --sequences x.fasta y.fasta \
  >   --output-sequences - > merged.fasta \
  >   --quiet
  ERROR: IDs must be unique within a sequence input file.
  
  The following entry is duplicated in 'y.fasta':
  
    'three'
  
  [2]

FASTA files without trailing newlines are supported.

  $ cat >y.fasta <<~~
  > >two
  > ATCG
  > >three
  > GCTA
  > ~~

  $ truncate -s -1 x.fasta
  $ truncate -s -1 y.fasta

  $ ${AUGUR} merge \
  >   --sequences x.fasta y.fasta \
  >   --output-sequences - > merged.fasta \
  >   --quiet

  $ cat merged.fasta
  >two
  ATCG
  >three
  GCTA
  >one
  ATCG
