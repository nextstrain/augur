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

The 'region name' column is query-able by backtick quoting.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query '(`region name` == "A")' \
  >  --output-strains filtered_strains.txt 2>/dev/null

  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
