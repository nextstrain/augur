Setup

  $ source "$TESTDIR"/_setup.sh

Using a query with a nonexistent column results in a specific error.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-pandas "invalid == 'value'" \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: Column 'invalid' does not exist in the metadata file. Ignoring it.
  ERROR: Query contains a column that does not exist in metadata.
  [2]

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-sqlite "invalid = 'value'" \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: Column 'invalid' does not exist in the metadata file. Ignoring it.
  ERROR: Query contains a column that does not exist in metadata.
  [2]


Using bad syntax in some queries results in meaningful errors, so they are exposed:

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-pandas "region >= 0.50" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Internal Pandas error when applying query:
  	'>=' not supported between instances of 'str' and 'float'
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]

FIXME: SQLite is not strongly typed, so this does not result in a syntax error:

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-sqlite "region >= 0.50" \
  >  --output-strains filtered_strains.txt
  0 strains were dropped during filtering
  12 strains passed all filters

However, other errors are not so helpful, so a link is provided for users to learn more about query syntax.

Unlike SQLite, Pandas does not understand '='.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-pandas "virus = 'zika'" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Internal Pandas error when applying query:
  	cannot assign without a target object
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]

Nonsensical queries behave the same for both Pandas and SQLite.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-pandas "some bad syntax" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Could not infer columns from the pandas query. If the query is valid, please specify columns using --query-columns.
  [2]

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query-sqlite "some bad syntax" \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: Column 'bad syntax' does not exist in the metadata file. Ignoring it.
  ERROR: Error when applying query. Ensure the syntax is valid per <https://www.sqlite.org/lang_expr.html>.
  [2]
