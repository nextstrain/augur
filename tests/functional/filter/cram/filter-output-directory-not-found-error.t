Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to output to a directory that does not exist.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --output-strains "directory-does-not-exist/filtered_strains.txt" > /dev/null
  ERROR: No such file or directory: 'directory-does-not-exist/filtered_strains.txt'
  [1]
