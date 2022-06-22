Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh


Test that a force-included strain is only output once.


Create some files for testing.

  $ cat >$TMP/metadata.tsv <<~~
  > strain	col
  > a	1
  > b	2
  > c	3
  > d	4
  > ~~
  $ cat >$TMP/sequences.fasta <<~~
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
  >   --metadata $TMP/metadata.tsv \
  >   --sequences $TMP/sequences.fasta \
  >   --subsample-max-sequences 4 \
  >   --include-where col=1 \
  >   --subsample-seed 0 \
  >   --output-metadata $TMP/metadata-filtered.tsv \
  >   --output-strains $TMP/strains-filtered.txt \
  >   --output-sequences $TMP/sequences-filtered.fasta \
  >   > /dev/null 2>&1
  $ cat $TMP/metadata-filtered.tsv | tail -n+2 | sort -k1
  a\t1 (esc)
  b\t2 (esc)
  c\t3 (esc)
  d\t4 (esc)
  $ cat $TMP/strains-filtered.txt | sort
  a
  b
  c
  d
  $ cat $TMP/sequences-filtered.fasta
  >a
  NNNN
  >b
  NNNN
  >c
  NNNN
  >d
  NNNN

Test all outputs with --include.

  $ cat >$TMP/include.txt <<~~
  > a
  > ~~
  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata.tsv \
  >   --sequences $TMP/sequences.fasta \
  >   --subsample-max-sequences 4 \
  >   --include $TMP/include.txt \
  >   --subsample-seed 0 \
  >   --output-metadata $TMP/metadata-filtered.tsv \
  >   --output-strains $TMP/strains-filtered.txt \
  >   --output-sequences $TMP/sequences-filtered.fasta \
  >   > /dev/null 2>&1
  $ cat $TMP/metadata-filtered.tsv | tail -n+2 | sort -k1
  a\t1 (esc)
  b\t2 (esc)
  c\t3 (esc)
  d\t4 (esc)
  $ cat $TMP/strains-filtered.txt | sort
  a
  b
  c
  d
  $ cat $TMP/sequences-filtered.fasta
  >a
  NNNN
  >b
  NNNN
  >c
  NNNN
  >d
  NNNN
