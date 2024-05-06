Setup

  $ source "$TESTDIR"/_setup.sh

  $ cat >metadata-no-date.tsv <<~~
  > strain	col1	col2	col3
  > SEQ1			b
  > SEQ2			b
  > SEQ3		c	d
  > SEQ4		c	d
  > ~~

An empty value in a --group-by column is still treated as a value for grouping.

I.e. the groups here are:
1. (None, None, b)
2. (None, c   , d)

  $ ${AUGUR} filter \
  >   --metadata metadata-no-date.tsv \
  >   --group-by col1 col2 col3 \
  >   --sequences-per-group 1 \
  >   --subsample-seed 0 \
  >   --output-log filtered-log.tsv \
  >   --output-strains filtered-strains.txt 2>/dev/null
  $ cat filtered-strains.txt
  SEQ1
  SEQ3
  $ tail -n+2 filtered-log.tsv | sort
  SEQ2\tsubsampling\t (esc)
  SEQ4\tsubsampling\t (esc)
