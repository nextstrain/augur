Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Error on duplicates in metadata.

Pandas engine
-------------

Within same chunk:

  $ cat >$TMP/metadata-duplicates.tsv <<~~
  > strain	date
  > a	2010-10-10
  > a	2010-10-10
  > b	2010-10-10
  > c	2010-10-10
  > d	2010-10-10
  > ~~
  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 10 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: Duplicate found in .* (re)
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

Separate chunks:

  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 1 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: Duplicate found in .* (re)
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >   --metadata $TMP/metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: Duplicate found in .* (re)
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]