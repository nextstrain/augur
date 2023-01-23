Setup

  $ source "$TESTDIR"/_setup.sh

Using a pandas query with a nonexistent column results in a specific error.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "invalid == 'value'" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Query contains a column that does not exist in metadata.
  [2]


Using pandas queries with bad syntax results in a generic errors.

This raises a ValueError internally (https://github.com/nextstrain/augur/issues/940):

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "invalid = 'value'" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Error when applying query. Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]

This raises a SyntaxError internally (https://github.com/nextstrain/augur/issues/941):

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "some bad syntax" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Error when applying query. Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]
