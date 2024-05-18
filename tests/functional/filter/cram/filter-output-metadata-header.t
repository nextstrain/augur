Setup

  $ source "$TESTDIR"/_setup.sh

Since Pandas's read_csv() and to_csv() are used with a double-quote character as
the default quotechar, any column names with that character may be altered.

Quoted columns containing the tab delimiter are left unchanged.

# FIXME: tsv-join has different behavior here. Test both?

  $ cat >metadata.tsv <<~~
  > strain	"col	1"
  > SEQ_1	a
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-metadata filtered_metadata.tsv 2>/dev/null

  $ head -n 1 filtered_metadata.tsv
  strain	"col	1"

Quoted columns without the tab delimiter are stripped of the quotes.

  $ cat >metadata.tsv <<~~
  > strain	"col1"
  > SEQ_1	a
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-metadata filtered_metadata.tsv 2>/dev/null

  $ head -n 1 filtered_metadata.tsv
  strain	col1

Any other columns with quotes are quoted, and pre-existing quotes are escsaped by doubling up.

  $ cat >metadata.tsv <<~~
  > strain	col"1	col2"
  > SEQ_1	a	b
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-metadata filtered_metadata.tsv 2>/dev/null

  $ head -n 1 filtered_metadata.tsv
  strain	"col""1"	"col2"""
