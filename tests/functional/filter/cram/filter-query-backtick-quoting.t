Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	region name
  > SEQ_1	A
  > SEQ_2	A
  > SEQ_3	B
  > SEQ_4	
  > ~~

The 'region name' column should be query-able by backtick quoting.
This does not currently work due to a bug.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query '(`region name` == "A")' \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: Could not infer columns from the pandas query. If the query is valid, please specify columns using --query-columns.
  [2]

  $ sort filtered_strains.txt
  sort: No such file or directory
  [2]
