Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	coverage	category
  > SEQ_1	0.94	A
  > SEQ_2	0.95	B
  > SEQ_3	0.96	C
  > SEQ_4		
  > ~~

The 'coverage' column should be query-able by numerical comparisons.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95" \
  >  --output-strains filtered_strains.txt > /dev/null

  $ sort filtered_strains.txt
  SEQ_2
  SEQ_3

The 'category' column will fail when used with a numerical comparison.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "category >= 0.95" \
  >  --output-strains filtered_strains.txt
  ERROR: Internal Pandas error when applying query:
  	'>=' not supported between instances of 'str' and 'float'
  Ensure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.
  [2]
