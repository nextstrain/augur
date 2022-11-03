Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Using a pandas query with a nonexistent column results in a specific error.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --query "invalid == 'value'" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: Query contains a column that does not exist in metadata.
  [1]


Using pandas queries with bad syntax results in a generic errors.

This raises a ValueError internally (https://github.com/nextstrain/augur/issues/940):

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --query "invalid = 'value'" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: Error when applying query. Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [1]

This raises a SyntaxError internally (https://github.com/nextstrain/augur/issues/941):

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --query "some bad syntax" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: Error when applying query. Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [1]
