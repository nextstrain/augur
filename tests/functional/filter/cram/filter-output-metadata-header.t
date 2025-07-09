Setup

  $ source "$TESTDIR"/_setup.sh

Quoting is unchanged regardless of placement.

  $ cat >metadata.tsv <<~~
  > strain	"col	1"
  > SEQ_1	a
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-metadata filtered_metadata.tsv 2>/dev/null

  $ head -n 1 filtered_metadata.tsv
  strain	"col	1"

  $ cat >metadata.tsv <<~~
  > strain	"col1"
  > SEQ_1	a
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-metadata filtered_metadata.tsv 2>/dev/null

  $ head -n 1 filtered_metadata.tsv
  strain	"col1"

  $ cat >metadata.tsv <<~~
  > strain	col"1	col2"
  > SEQ_1	a	b
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-metadata filtered_metadata.tsv 2>/dev/null

  $ head -n 1 filtered_metadata.tsv
  strain	col"1	col2"
