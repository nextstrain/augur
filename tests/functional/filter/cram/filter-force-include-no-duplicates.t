Setup

  $ source "$TESTDIR"/_setup.sh


Test that a force-included strain is only output once.


Create some files for testing.

  $ cat >metadata.tsv <<~~
  > strain	col
  > a	1
  > b	2
  > c	3
  > d	4
  > ~~
  $ cat >sequences.fasta <<~~
  > >a
  > NNNN
  > >b
  > NNNN
  > >c
  > NNNN
  > >d
  > NNNN
  > ~~

Test all outputs with --include-where.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --subsample-max-sequences 4 \
  >   --include-where col=1 \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv \
  >   --output-strains strains-filtered.txt \
  >   --output-sequences sequences-filtered.fasta \
  >   2>/dev/null
  $ cat metadata-filtered.tsv | tail -n+2 | sort -k1
  a\t1 (esc)
  b\t2 (esc)
  c\t3 (esc)
  d\t4 (esc)
  $ cat strains-filtered.txt | sort
  a
  b
  c
  d
  $ cat sequences-filtered.fasta
  >a
  NNNN
  >b
  NNNN
  >c
  NNNN
  >d
  NNNN

Test all outputs with --include.

  $ cat >include.txt <<~~
  > a
  > ~~
  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --subsample-max-sequences 4 \
  >   --include include.txt \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv \
  >   --output-strains strains-filtered.txt \
  >   --output-sequences sequences-filtered.fasta \
  >   2>/dev/null
  $ cat metadata-filtered.tsv | tail -n+2 | sort -k1
  a\t1 (esc)
  b\t2 (esc)
  c\t3 (esc)
  d\t4 (esc)
  $ cat strains-filtered.txt | sort
  a
  b
  c
  d
  $ cat sequences-filtered.fasta
  >a
  NNNN
  >b
  NNNN
  >c
  NNNN
  >d
  NNNN
