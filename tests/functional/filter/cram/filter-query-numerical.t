Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	coverage
  > SEQ_1	0.94
  > SEQ_2	0.95
  > SEQ_3	0.96
  > SEQ_4	
  > ~~

Ideally, the 'coverage' column should be query-able by numerical comparisons.
This does not currently work since the empty string is causing that column to be
parsed as a non-numerical type.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95" \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Internal Pandas error when applying query:
  	'>=' not supported between instances of 'str' and 'float'
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]
