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

The 'coverage' column should be query-able by numerical comparisons.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95" \
  >  --output-strains filtered_strains.txt > /dev/null

  $ sort filtered_strains.txt
  SEQ_2
  SEQ_3
