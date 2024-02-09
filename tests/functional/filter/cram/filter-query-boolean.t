Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	True
  > SEQ_2	True
  > SEQ_3	False
  > ~~

Ideally, the column should be query-able by boolean comparisons.
This does not currently work because all dtypes are strings.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column == True" \
  >  --output-strains filtered_strains.txt
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  3 strains were dropped during filtering
  \t3 were filtered out by the query: "column == True" (esc)
  [2]

  $ sort filtered_strains.txt
