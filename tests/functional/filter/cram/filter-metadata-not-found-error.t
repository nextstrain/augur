Setup

  $ source "$TESTDIR"/_setup.sh

Try to filter on an metadata file that does not exist.

  $ ${AUGUR} filter \
  >  --metadata file-does-not-exist.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: No such file or directory: 'file-does-not-exist.tsv'
  [2]
