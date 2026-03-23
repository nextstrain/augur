SETUP

  $ source "$TESTDIR"/_setup.sh


SINGLE METADATA INPUT

A single metadata input is supported, e.g. to annotate rows with a source
column or to take advantage of compression support.

  $ cat >x.tsv <<~~
  > strain	a	b	c
  > one	X1a	X1b	X1c
  > two	X2a	X2b	X2c
  > ~~

Pass-through with source column.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv \
  >   --source-columns '__source_metadata_{NAME}' \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    __source_metadata_X
  one     X1a  X1b  X1c                    1
  two     X2a  X2b  X2c                    1

Pass-through without source column.

  $ ${AUGUR} merge \
  >   --metadata X=x.tsv \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c
  one     X1a  X1b  X1c
  two     X2a  X2b  X2c

Supports Augur's standard accepted compression formats.

  $ xz < x.tsv > x.tsv.xz
  $ ${AUGUR} merge \
  >   --metadata X=x.tsv.xz \
  >   --source-columns '__source_metadata_{NAME}' \
  >   --output-metadata - --quiet | csv2tsv --csv-delim $'\t' | tsv-pretty
  strain  a    b    c    __source_metadata_X
  one     X1a  X1b  X1c                    1
  two     X2a  X2b  X2c                    1


SINGLE SEQUENCE INPUT

A single sequence input is supported, e.g. to take advantage of compression
support.

  $ cat >x.fasta <<~~
  > >seq1
  > ATCG
  > >seq2
  > GCTA
  > ~~

  $ ${AUGUR} merge \
  >   --sequences x.fasta \
  >   --output-sequences merged.fasta
  Validating 'x.fasta'…
  Merging sequences and writing to 'merged.fasta'…

  $ cat merged.fasta
  >seq1
  ATCG
  >seq2
  GCTA
