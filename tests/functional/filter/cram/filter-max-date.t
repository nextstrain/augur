Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	2020-03-XX
  > SEQ_2	2020-03-01
  > SEQ_3	2020-03-02
  > ~~

Test that --max-date is inclusive.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --max-date 2020-03-01 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
