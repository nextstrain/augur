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

Create another metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	metric1	metric2
  > SEQ1	4	5
  > SEQ2	5	9
  > SEQ3	6	10
  > ~~

Use a Pandas query to filter by a numerical value.
This relies on having proper data types associated with the columns. If < is
comparing strings, it's likely that SEQ3 will be dropped or errors arise.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "metric1 > 4 & metric1 < metric2" \
  >  --output-strains filtered_strains.txt
  1 strains were dropped during filtering
  \t1 were filtered out by the query: "metric1 > 4 & metric1 < metric2" (esc)
  2 strains passed all filters

  $ sort filtered_strains.txt
  SEQ2
  SEQ3
