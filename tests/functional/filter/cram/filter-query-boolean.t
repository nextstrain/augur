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

Note that the string value is case-insensitive.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	True
  > SEQ_2	trUe
  > SEQ_3	FALSE
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column == True" \
  >  --output-strains filtered_strains.txt
  1 strain was dropped during filtering
  	1 was filtered out by the query: "column == True"
  2 strains passed all filters

Note that 1/0 can also be compared to boolean literals.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	1
  > SEQ_2	1
  > SEQ_3	0
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column == True" \
  >  --output-strains filtered_strains.txt
  1 strain was dropped during filtering
  	1 was filtered out by the query: "column == True"
  2 strains passed all filters

Empty values are ignored.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	True
  > SEQ_2	False
  > SEQ_3	
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column == True" \
  >  --output-strains filtered_strains.txt
  2 strains were dropped during filtering
  	2 were filtered out by the query: "column == True"
  1 strain passed all filters

  $ sort filtered_strains.txt
  SEQ_1
