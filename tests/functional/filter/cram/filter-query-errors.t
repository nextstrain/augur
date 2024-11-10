Setup

  $ source "$TESTDIR"/_setup.sh

Using a pandas query with a nonexistent column results in a specific error.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "invalid == 'value'" \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: Column 'invalid' does not exist in the metadata file. This may cause subsequent errors.
  ERROR: Query contains a column that does not exist in metadata: name 'invalid' is not defined
  [2]


Using pandas queries with bad syntax results in meaningful errors.

Some error messages from Pandas may be useful, so they are exposed:

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "region >= 0.50" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Internal Pandas error when applying query:
  	'>=' not supported between instances of 'str' and 'float'
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]

However, other Pandas errors are not so helpful, so a link is provided for users to learn more about query syntax.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "country = 'value'" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Internal Pandas error when applying query:
  	cannot assign without a target object
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "some bad syntax" \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: Could not infer columns from the pandas query. Reading all metadata columns,
  which may impact execution time. If the query is valid, please open a new issue:
  
      <https://github.com/nextstrain/augur/issues/new/choose>
  
  and add the query to the description:
  
      some bad syntax
  ERROR: Internal Pandas error when applying query:
  	invalid syntax (<unknown>, line 1)
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]
