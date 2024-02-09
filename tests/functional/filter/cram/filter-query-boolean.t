Setup

  $ source "$TESTDIR"/_setup.sh

A column with True and False values is query-able by boolean comparisons.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	True
  > SEQ_2	True
  > SEQ_3	False
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column == True" \
  >  --output-strains filtered_strains.txt
  1 strain was dropped during filtering
  	1 was filtered out by the query: "column == True"
  2 strains passed all filters

  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
