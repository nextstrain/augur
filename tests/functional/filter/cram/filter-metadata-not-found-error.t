Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to filter on an metadata file that does not exist.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --metadata file-does-not-exist.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: No such file or directory: 'file-does-not-exist.tsv'
  [2]

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata file-does-not-exist.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: No such file or directory: 'file-does-not-exist.tsv'
  [2]
