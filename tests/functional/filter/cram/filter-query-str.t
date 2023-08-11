Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	value1
  > SEQ_2	value2
  > SEQ_3	value3
  > ~~

'column' should be query-able using the `.str` accessor.
Currently, this does not work properly due to a bug.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column.str.startswith('value')" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Internal Pandas error when applying query:
  	'NoneType' object has no attribute 'str'
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]
