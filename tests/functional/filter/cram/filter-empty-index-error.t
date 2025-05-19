Setup

  $ source "$TESTDIR"/_setup.sh

Error on empty indexes in metadata.

  $ cat >metadata-empty-indexes.tsv <<~~
  > strain	date
  > 	2010-10-10
  > 	2010-10-10
  > b	2010-10-10
  > c	2010-10-10
  > d	2010-10-10
  > ~~
  $ ${AUGUR} filter \
  >   --metadata metadata-empty-indexes.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 10 \
  >   --output-metadata metadata-filtered.tsv > /dev/null
  ERROR: Found rows with empty values in id column 'strain' in .* (re)
  Please remove the rows with empty ids or use a different id column via --metadata-id-columns.
  [2]
  $ cat metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]
