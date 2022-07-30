Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Error on missing group-by columns.

Pandas engine
-------------

  $ cat >$TMP/metadata-no-date.tsv <<~~
  > strain	col
  > SEQ1	a
  > SEQ2	b
  > ~~

  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata-no-date.tsv \
  >   --group-by year \
  >   --sequences-per-group 1 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: The specified group-by categories (['year']) were not found. Note that using 'year' or 'year month' requires a column called 'date'.
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata-no-date.tsv \
  >   --group-by invalid \
  >   --sequences-per-group 1 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: The specified group-by categories (['invalid']) were not found.
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

SQLite engine
-------------

  $ cat >$TMP/metadata-no-date.tsv <<~~
  > strain	col
  > SEQ1	a
  > SEQ2	b
  > ~~

  $ ${AUGUR} filter --engine sqlite \
  >   --metadata $TMP/metadata-no-date.tsv \
  >   --group-by year \
  >   --sequences-per-group 1 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: The specified group-by categories (['year']) were not found. Note that using 'year' or 'year month' requires a column called 'date'.
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

  $ ${AUGUR} filter --engine sqlite \
  >   --metadata $TMP/metadata-no-date.tsv \
  >   --group-by invalid \
  >   --sequences-per-group 1 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: The specified group-by categories (['invalid']) were not found.
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]
