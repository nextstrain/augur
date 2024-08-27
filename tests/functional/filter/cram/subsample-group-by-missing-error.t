Setup

  $ source "$TESTDIR"/_setup.sh

Error on missing group-by columns.

  $ cat >metadata-no-date.tsv <<~~
  > strain	col
  > SEQ1	a
  > SEQ2	b
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata-no-date.tsv \
  >   --group-by year \
  >   --sequences-per-group 1 \
  >   --output-metadata metadata-filtered.tsv > /dev/null
  WARNING: Column 'date' does not exist in the metadata file. This may cause subsequent errors.
  ERROR: The specified group-by categories (['year']) were not found. Note that using any of ['month', 'week', 'year'] requires a column called 'date'.
  [2]
  $ cat metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

  $ ${AUGUR} filter \
  >   --metadata metadata-no-date.tsv \
  >   --group-by invalid \
  >   --sequences-per-group 1 \
  >   --output-metadata metadata-filtered.tsv > /dev/null
  WARNING: Column 'invalid' does not exist in the metadata file. This may cause subsequent errors.
  ERROR: The specified group-by categories (['invalid']) were not found.
  [2]
  $ cat metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]
